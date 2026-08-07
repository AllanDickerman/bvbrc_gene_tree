package Tree_Builder;
use Sequence_Alignment;
use IPC::Run 'run';
use File::Copy ('copy', 'move');
use File::Path qw(make_path remove_tree);
use File::Temp;
use File::Basename;
use Cwd ('abs_path', 'getcwd');
use strict;
use warnings;
our $debug = 0;

sub set_debug {
    $debug = shift;
}

sub new {
    my ($class, $alignment) = @_;
    if ($debug) {
        print(STDERR "in Tree_Builder\n");
        print(STDERR " class = $class\n");
        print(STDERR " alignment = $alignment\n");
    }
    my $self = {};
    bless $self, $class;
    $self->{_alignment_file} = $alignment;
    my $output_base = $alignment;
    $output_base =~ s/.*\///; # remove path if it exists
    $output_base =~ s/^(.*)\..*/$1/; # remove extension if it exists
    print STDERR "output_base = $output_base\n";
    $self->{_output_base} = $output_base;

    $self->{_parallel} = 2;
    $self->{_model} = 'GTR';
    $self->{_program} = 'fasttree'; # or raxml or phyml
    $self->{_analysis_steps} = ();
    $self->{_original_dir} = getcwd(); 
    return $self;
}

sub set_output_base {
    my ($self, $output_base) = @_;
    $self->{_output_base} = $output_base;
}

sub set_output_dir {
    my ($self, $output_dir) = @_;
    unless ( -d $output_dir ) {
        make_path $output_dir;
    }
    $self->{_output_dir} = $output_dir;
}

sub set_alphabet {
    my ($self, $alphabet) = @_;
    $self->{_alphabet} = $alphabet;
}

sub set_model {
    my ($self, $model) = @_;
    $self->{_model} = $model;
}

sub set_program { # raxml, phyml, or fasttree
    my ($self, $program) = @_;
    $self->{_program} = lc($program);
}

sub autodetect_alphabet {
    my $alignment_file = shift;
    my $dna_count = 0;
    my $protein_count = 0;
    open F, $alignment_file or warn "cannot open $alignment_file";
    while (<F>) {
        next if /^>/;
        s/^\S+\s+(\S.*)/$1/; # remove taxon name in case of phylip format (for phyml)
        $dna_count += tr/acgtACGT/acgtACGT/;
        $protein_count += tr/defhiklmnopqrsuvwyDEFHIKLMNOPQRSUVWY/defhiklmnopqrsuvwyDEFHIKLMNOPQRSUVWY/;
    }
    print STDERR "autodetect_alphabet counted $dna_count DNA letters, $protein_count protein letters.\n" if $debug;
    return ('DNA', 'protein')[$protein_count > $dna_count];
}

sub get_analysis_steps {
    my $self = shift;
    return $self->{_analysis_steps};
}

sub get_analysis_name {
    my $self = shift;
    my $retval;
    for my $step (@{$self->{_analysis_steps}}) {
        if ($retval) {
            $retval .= "\n";
        }
        $retval .= $step->{name};
    }
    return $retval;
}

sub get_analysis_commandline {
    my $self = shift;
    my $retval;
    for my $step (@{$self->{_analysis_steps}}) {
        if ($retval) {
            $retval .= "\n";
        }
        $retval .= $step->{command_line};
    }
    return $retval;
}

sub get_analysis_stderrout {
    my $self = shift;
    my $retval;
    for my $step (@{$self->{_analysis_steps}}) {
        if ($retval) {
            $retval .= "\n";
        }
        $retval .= $step->{stdout};
        $retval .= $step->{stderr};
    }
    return $retval;
}

sub get_analysis_names {
    my $self = shift;
    my $retval;
    for my $step (@{$self->{_analysis_steps}}) {
        $retval .= $step->{name};
    }
    return $retval;
}

sub add_analysis_step {
    my ($self, $descriptor, $command_line) = @_;
    my %step;
    $step{name} = $descriptor;
    $step{command_line} = $command_line;
    $step{start_time} = time();
    $step{comments} = [];
    print STDERR "$descriptor\n$command_line\n" if $debug;
    push @{$self->{_analysis_steps}}, \%step;
}

sub add_analysis_out_err {
    my ($self, $stdout, $stderr) = @_;
    my $index = $#{$self->{_analysis_steps}};
    my $step = $self->{_analysis_steps}[$index];
    $step->{stdout} = $stdout;
    $step->{stderr} = $stderr;
    $step->{end_time} = time();
}

sub get_step_info {
    my $self = shift;
    my $index = $#{$self->{_analysis_steps}};
    my $step = $self->{_analysis_steps}[$index];
    return $step;
}
    
sub add_analysis_comment {
    my ($self, $comment) = @_;
    my $index = $#{$self->{_analysis_steps}};
    my $step = $self->{_analysis_steps}[$index];
    push @{$step->{comments}}, $comment;
}
sub add_analysis_tree {
    my ($self, $tree_file) = @_;
    my $index = $#{$self->{_analysis_steps}};
    my $step = $self->{_analysis_steps}[$index];
    $step->{tree_file} = $tree_file;
}
sub add_analysis_log {
    my ($self, $log_file) = @_;
    my $index = $#{$self->{_analysis_steps}};
    my $step = $self->{_analysis_steps}[$index];
    $step->{log_file} = $log_file;
}

sub get_num_analysis_steps {
    my $self = shift;
    return scalar @{$self->{_analysis_steps}};
} 

sub get_analysis_step {
    my ($self, $index) = @_;
    if ($index > $#{$self->{_analysis_steps}}) {
        print STDERR "Problem: index passed to get_analysis_step, $index, is higher than available steps.\n";
        return;
    }
    return $self->{_analysis_steps}[$index];
} 

sub get_log_file {
    my $self = shift;
    my $index = $#{$self->{_analysis_steps}};
    my $retval = $self->{_analysis_steps}[$index]->{log_file};
    return $retval;
} 

sub clear_analysis_history {
    my $self = shift;
    $self->{_analysis_steps} = []
}

sub set_bootstrap_reps {
    my ($self, $boostrap) = @_;
    $self->{_bootstrap} = $boostrap;
}

sub build_tree {
    my ($self, $program) = @_;
    $program = "" unless $program;
    print "build_tree($program)\n";
    if ($program) {
        $self->set_program($program);
    }
    if (! $self->{_alphabet} ) {
        $self->{_alphabet} = autodetect_alphabet($self->{_alignment_file});
    }

    $self->{_tmpdir} = File::Temp->newdir( "/tmp/Tree_Builder_XXXXX", CLEANUP => !$debug );
    system("chmod", "755", $self->{_tmpdir});
    print STDERR "created temp dir: $self->{_tmpdir}, cleanup = ", !$debug, "\n";
    my $alignment_abs_path = abs_path($self->{_alignment_file});
    print  STDERR "Alignment abs path = $alignment_abs_path\n" if $debug;
    chdir($self->{_tmpdir}); # do all work in temporary directory
    my $rv = symlink($alignment_abs_path, $self->{_alignment_file});
    print("symlink($alignment_abs_path, $self->{_alignment_file})\n");
    unless ($rv) {
        print("return val of symlink = $rv\n");
        chdir($self->{_original_dir});
        die "symlink failed: $!" unless $rv;
    }
    system("ls -l $self->{_tmpdir}");
    my $treefile;
    system("pwd");
    system("ls -l ");
    print("program = $self->{_program}\n");
    if ($self->{_program} eq 'fasttree') {
        $treefile = $self->build_fasttree();
    }
    elsif ($self->{_program} eq 'raxml') {
       $treefile =  $self->build_raxml_tree();
    }
    elsif ($self->{_program} eq 'phyml') {
        $treefile = $self->build_phyml_tree();
    }
    else {
        chdir($self->{_original_dir}); # go back
        die("cannot interpret program: $self->{_program}");
    }
    chdir($self->{_original_dir}); # go back
    return $treefile;
}

sub build_raxml_tree {
    my $self = shift;
    
    $self->clear_analysis_history();
    #$self->{_program} = 'raxml';
    print STDERR "In build_raxml_tree\n" if $debug;

    my $model = uc($self->{_model}); 
    if ($self->{_alphabet} eq 'DNA') {
        $model = 'GTRGAMMA'
    }
    else {
        $model = 'LG' if $model !~ /DAYHOFF|DCMUT|JTT|MTREV|WAG|RTREV|CPREV|VT|BLOSUM62|MTMAM|LG|MTART|MTZOA|PMB|HIVB|HIVW|JTTDCMUT|FLU|STMTREV|DUMMY|DUMMY2|AUTO|LG4M|LG4X|PROT_FILE/i;
        $model = "PROTCAT". $model;
    }
    #$self->{_model} = $model;

    my @cmd = ("raxmlHPC-PTHREADS-SSE3");
    push @cmd, ("-T", $self->{_parallel});
    push @cmd, ("-p", "12345");
    
    push @cmd, ("-x", "12345");
    push @cmd, ("-N", "100");
    push @cmd, ("-m", $model);
    push @cmd, ("-s", $self->{_alignment_file});
    push @cmd, ("-n", $self->{_output_base});
    push @cmd, ("-f", "a"); # maximum-likelihood tree with rapid boostrap support values
    
    my $analysis_descriptor = "Build maximum-likelihood tree using RAxML";
    $self->add_analysis_step($analysis_descriptor, join(" ", @cmd));
    print STDOUT "run command: ". join(" ", @cmd)."\n" if $debug;
    my ($out, $err) = run_cmd(\@cmd);
    $self->add_analysis_out_err($out, $err);
    my $tree_name = $self->{_output_base} . "_raxml.nwk";
    rename("RAxML_bipartitions." . $self->{_output_base}, $tree_name);
    $self->add_analysis_tree($tree_name);
    move("RAxML_info." . $self->{_output_base}, $self->{_original_dir});
    move($tree_name, $self->{_original_dir});
    return $tree_name;
}

sub build_phyml_tree {
    my $self = shift;
    my $datatype = 'aa';
    my $model = $self->{_model};
    print("build_phyml_tree(), model=$model)\n");
    if ($self->{_alphabet} =~ /DNA/i) {
        $datatype = 'nt';
        $model = 'GTR' unless $model =~ /HKY85|JC69|K80|F81|F84|TN93|GTR/;
    }
    else {
        $model = 'LG' unless $model =~ /WAG|JTT|MtREV|Dayhoff|DCMut|RtREV|CpREV|VT|AB|Blosum62|MtMam|MtArt|HIVw|HIVb/;
    }
    print("build_phyml_tree(), model=$model)\n");

    my $analysis_descriptor = "Build maximum-likelihood tree using PhyML";
    my @cmd = ("phyml");
    push @cmd, ("-i", $self->{_alignment_file});
    push @cmd, ("-d", $datatype);
    push @cmd, ("-m", $model);
    my $comment;
    if ($self->{_bootstrap}) {
        push @cmd, ("-b", $self->{_bootstrap}); # normal bootstrap replicates
        $comment = "With $self->{_bootstrap} bootstrap replicates.\n";
        print STDERR $comment if $debug;
    }
    else {
        push @cmd, ("-b", '-5'); # -5 gives approximate Bayes branch supports.
        # -2 give Chi2-based parametric branch supports
        # -4 gives Shimodaira-Hasegawa (SH) support values
        $comment = "Defaults to approximate Bayes branch supports.\n";
        print STDERR $comment if $debug;
    }
    if ($ENV{P3_ALLOCATED_CPU}) {
        # set PHYMLCPUS environment variable
        print STDOUT "set PHYMLCPUS to $ENV{P3_ALLOCATED_CPU}\n";
        $ENV{PHYMLCPUS} = $ENV{P3_ALLOCATED_CPU};
    }
   
    
    $self->add_analysis_step($analysis_descriptor, join(" ", @cmd));
    $self->add_analysis_comment($comment);
    print STDOUT "run phyml command: ". join(" ", @cmd)."\n" if $debug;
    my ($out, $err) = run_cmd(\@cmd);
    $self->add_analysis_out_err($out, $err);
    print STDOUT "phmyl stdout:\n$out\n";

    my $treeFile = $self->{_output_base} . "_phyml_tree.nwk";
    move($self->{_alignment_file}."_phyml_tree.txt", "$self->{_original_dir}/$treeFile");# copy final tree to original working directory
    $self->add_analysis_tree($treeFile);
    my $logFile = $self->{_output_base} . "_phyml_log.txt";
    move($self->{_alignment_file}."_phyml_stats.txt", "$self->{_original_dir}/$logFile");
    $self->add_analysis_log($logFile);
    return $treeFile;
} 
    
sub build_fasttree {
    my $self = shift;
    my $treeFile = $self->{_output_base} . "_fasttree.nwk";
    my $logFile = $self->{_output_base} . "_fasttree_log.txt";
    my @cmd = ("FastTree", "-out", $treeFile, "-log", $logFile);
    if ($self->{_alphabet} =~ /DNA/i) {
        push @cmd, "-nt", "-gtr", "-gamma";
    }
    elsif ($self->{_model} =~ /lg|wag/i) {
        push @cmd, "-" . lc($self->{_model});
    } # else defaults to JTT for proteins

    push @cmd, $self->{_alignment_file};
    my $analysis_descriptor = "Build maximum-likelihood tree using FastTree"; 
    my $tree_file_name;
    $self->add_analysis_step($analysis_descriptor, join(" ", @cmd));
    print STDOUT "run command: ". join(" ", @cmd)."\n" if $debug;
    my ($out, $err) = run_cmd(\@cmd);
    $self->add_analysis_out_err($out, $err);
    $self->add_analysis_tree($treeFile);
    $self->add_analysis_log($logFile);

    move($treeFile, "$self->{_original_dir}/$treeFile");
    move($logFile, "$self->{_original_dir}/$logFile");
    return $treeFile;
}

sub run_cmd {
    my ($cmd) = @_;
    my ($out, $err);
    run($cmd, '>', \$out, '2>', \$err)
        or die "Error running cmd=@$cmd, stdout:\n$out\nstderr:\n$err\n";
    # print STDERR "STDOUT:\n$out\n";
    # print STDERR "STDERR:\n$err\n";
    return ($out, $err);
}


1
