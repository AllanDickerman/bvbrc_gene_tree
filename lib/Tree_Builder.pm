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
    my ($class, $fasta_alignment, $alphabet, $threads) = @_;
    if ($debug) {
        print(STDERR "in Tree_Builder\n");
        print(STDERR " class = $class\n");
        print(STDERR " fasta_alignment = $fasta_alignment\n");
    }
    my $self = {};
    bless $self, $class;
    my $alignment_abs_path = File::Spec->rel2abs( $fasta_alignment );
    $self->{_alignment_file} = $fasta_alignment;
    $self->{_alphabet} = $alphabet;
    $self->{_parallel} = $threads ? $threads : 1;
    my ($fasta_name, $fasta_dir, $fasta_extension) =  fileparse($alignment_abs_path, qr/\.[^.]*/);
    $self->{_output_base} = $fasta_name;
    $self->{_alignment_name} = $fasta_name . $fasta_extension;
    print STDERR "fileparse on $fasta_alignment yields $fasta_dir, $fasta_name, $fasta_extension\n" if $debug;
    $alphabet = autodetect_alphabet($fasta_alignment) unless $alphabet;
    $self->{_model} = ('GTR', 'LG')[$alphabet eq 'protein']; #default to GTR for DNA and LG for proteins, can change
    $self->{_analysis_steps} = ();
    $self->{_tmpdir} = File::Temp->newdir( "/tmp/Tree_Builder_XXXXX", CLEANUP => !$debug );
    system("chmod", "755", $self->{_tmpdir});
    print STDERR "created temp dir: $self->{_tmpdir}, cleanup = ", !$debug, "\n";
    $self->{_output_dir} = $self->{_original_wd} = getcwd(); 
    chdir($self->{_tmpdir}); # do all work in temporary directory
    symlink($alignment_abs_path, $self->{_alignment_name});
    return $self;
}

DESTROY {
    #destructor for Tree_Builder object
    my $self = shift;
    print STDERR "In Tree_Builder::DESTROY\n" if $debug;
    #chdir($self->{_original_wd}); # go back to original directory
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

sub set_model {
    my ($self, $model) = @_;
    $self->{_model} = $model;
}

sub autodetect_alphabet {
    my $alignment_file = shift;
    my $dna_count = 0;
    my $protein_count = 0;
    open F, $alignment_file or warn "cannot open $alignment_file";
    while (<F>) {
        next if /^>/;
        $dna_count += tr/acgtACGT/acgtACGT/;
        $protein_count += tr/defhiklmnopqrsuvwyDEFHIKLMNOPQRSUVWY/defhiklmnopqrsuvwyDEFHIKLMNOPQRSUVWY/;
    }
    print STDERR "autodetect_alphabet counted $dna_count DNA letters, $protein_count protein letters.\n" if $debug;
    return ('DNA', 'protein')[$protein_count > $dna_count];
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

sub build_raxml_tree {
    my ($self, $bootstrap) = @_;
    chdir($self->{_tmpdir}); # do all work in temporary directory
    my $output_base = $self->{_output_base};
    my $logFile = "${output_base}_raxml_log.txt";
    $self->clear_analysis_history();
    $self->{_program} = 'raxml';
    print STDERR "In build_raxml_tree\n" if $debug;

    if ($debug) {
        for my $key (sort keys %$self) {
            print STDERR "tb:$key = \t$self->{$key}\n";
        }
        print "\n";
    }
    
    #my $cwd = getcwd();
    my $model = uc($self->{_model}); 
    if ($self->{_alphabet} eq 'DNA') {
        $model = 'GTRGAMMA'
    }
    else {
        $model = 'LG' if $model !~ /DAYHOFF|DCMUT|JTT|MTREV|WAG|RTREV|CPREV|VT|BLOSUM62|MTMAM|LG|MTART|MTZOA|PMB|HIVB|HIVW|JTTDCMUT|FLU|STMTREV|DUMMY|DUMMY2|AUTO|LG4M|LG4X|PROT_FILE/i;
        $model = "PROTCAT". $model;
    }
    #$self->{_model} = $model;

    my $analysis_descriptor = "Build maximum-likelihood tree using RAxML";
    my @cmd = ("raxmlHPC-PTHREADS-SSE3");
    push @cmd, ("-T", $self->{_parallel});
    push @cmd, ("-p", "12345");
    push @cmd, ("-m", $model);
    push @cmd, ("-s", $self->{_alignment_file});
    push @cmd, ("-n", $output_base);
    my $support_trees;
    my $comment;
    my $tree_with_support_name = "${output_base}_raxml";
    if ($bootstrap) {
        push @cmd, ("-f",  "d"); # just do ML search (generate bootstrap trees separately later)
        $tree_with_support_name .= "_bootstrap_tree.nwk";
    }
    else {
        push @cmd, ("-f", "D"); # generate RELL replicate trees
        $support_trees = "RAxML_rellBootstrap." . $output_base;
        $comment = "RELL replicates";
        $tree_with_support_name .= "_rell_tree.nwk";
    }
    
    $self->add_analysis_step($analysis_descriptor, join(" ", @cmd));
    $self->add_analysis_comment($comment) if $comment;
    print STDOUT "run command: ". join(" ", @cmd)."\n" if $debug;
    my ($out, $err) = run_cmd(\@cmd);
    $self->add_analysis_out_err($out, $err);
    $self->add_analysis_tree("RAxML_bestTree.$output_base");
    move("RAxML_info.$output_base", $logFile);
    $self->add_analysis_log($logFile);

    if ($bootstrap) {
        if ($self->{_alphabet} eq 'DNA') {
            $model = 'GTRCAT'
        }
        my $bootstrap_name = "Tree_Builder_bootstrap";
        $analysis_descriptor = "Generate bootstrap trees";
        $support_trees = "RAxML_bootstrap." . $bootstrap_name;
        @cmd = ("raxmlHPC-PTHREADS-SSE3");
        push @cmd, ("-T", $self->{_parallel});
        push @cmd, ("-m", $model);
        push @cmd, ("-s", $self->{_alignment_file});
        push @cmd, ("-n", $bootstrap_name);
        push @cmd, ("-b",  "12345", "-#", $bootstrap, "-p", "12345");
        $self->add_analysis_step($analysis_descriptor, join(" ", @cmd));
        print STDOUT "run command: ". join(" ", @cmd)."\n" if $debug;
        ($out, $err) = run_cmd(\@cmd);
        $self->add_analysis_out_err($out, $err);
        $self->add_analysis_tree($support_trees);
        $self->add_analysis_comment("Number of bootstrap replicates = $bootstrap.");
        open CAT, ">>$logFile";
        open IN, "RAxML_info.$output_base";
        print CAT <IN>;
        close CAT;
        $self->add_analysis_log($logFile);
    }
    # now map replicate support numbers onto ML tree (from either bootstrap or RELL output)
    $analysis_descriptor = "Map replicate support numbers onto ML tree";
    @cmd = ("raxmlHPC-PTHREADS-SSE3", "-f", "b", "-m", "GTRCAT"); #to map RELL bootstrap support onto ML tree
    push @cmd, ("-t", "RAxML_bestTree.".$output_base);
    push @cmd, ("-z", $support_trees);
    push @cmd, ("-n", $tree_with_support_name);
    push @cmd, ("-T", $self->{_parallel}, "-p", "12345");
    $self->add_analysis_step($analysis_descriptor, join(" ", @cmd));
    print STDOUT "run command: ". join(" ", @cmd)."\n" if $debug;
    ($out, $err) = run_cmd(\@cmd);
    $self->add_analysis_out_err($out, $err);
    move("RAxML_bipartitions.$tree_with_support_name", "$self->{_output_dir}/$tree_with_support_name");
    $self->add_analysis_tree($tree_with_support_name);
    open CAT, ">>$logFile";
    open IN, "RAxML_info.$tree_with_support_name";
    print CAT <IN>;
    close CAT;
    $self->add_analysis_log($logFile);
    move($logFile, "$self->{_output_dir}/$logFile");
    chdir($self->{_original_wd}); 
    return $tree_with_support_name;
}

sub build_phyml_tree {
    my ($self, $bootstrap) = @_;
    chdir($self->{_tmpdir}); 
    my $output_base = $self->{_output_base};
    unless ($self->{_phylip_file}) {
        Sequence_Alignment::set_debug(1) if $debug;
        my $alignment = new Sequence_Alignment($self->{_alignment_file});
        $self->{_phylip_file} = $self->{_alignment_file};
        $self->{_phylip_file} =~ s/\..{2,5}$//;
        $self->{_phylip_file} .= ".phy";
        $alignment->write_phylip($self->{_phylip_file});
    }
    my $datatype = 'aa';
    my $model = $self->{_model};
    if ($self->{_alphabet} =~ /DNA/i) {
        $datatype = 'nt';
        $model = 'GTR' if $model !~ /HKY85|JC69|K80|F81|F84|TN93|GTR/;
    }
    else {
        $model = 'LG' if $model !~ /WAG|JTT|MtREV|Dayhoff|DCMut|RtREV|CpREV|VT|AB|Blosum62|MtMam|MtArt|HIVw|HIVb/;
    }

    my $treeFile = $output_base."_phyml_tree.nwk";
    my $analysis_descriptor = "Build maximum-likelihood tree using PhyML";
    my @cmd = ("phyml");
    push @cmd, ("-i", $self->{_phylip_file});
    push @cmd, ("-d", $datatype);
    push @cmd, ("-m", $model);
    my $comment;
    if ($bootstrap) {
        push @cmd, ("-b", $bootstrap); # normal bootstrap replicates
        $comment = "With $bootstrap bootstrap replicates.";
        $treeFile = $output_base."_phyml_bootstrap_tree.nwk";
    }
    else {
        push @cmd, ("-b", '-3'); # -1 gives approximate Likelihood Ratio Tests
        # -2 give Chi2-based parametric branch supports
        # -3 gives Shimodaira-Hasegawa (SH) support values
        $comment = "Support values from the approximate likelihood ratio test.";
        $treeFile = $output_base."_phyml_alr_tree.nwk";
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
    move($self->{_phylip_file}."_phyml_tree.txt", "$self->{_output_dir}/$treeFile");# copy final tree to original working directory
    $self->add_analysis_tree($treeFile);
    my $logFile = $output_base."_phyml_log.txt";
    move($self->{_phylip_file}."_phyml_stats.txt", "$self->{_output_dir}/$logFile");
    $self->add_analysis_log($logFile);
    chdir($self->{_original_wd}); 
    return $treeFile;
} 
    
sub build_fasttree {
    my ($self, $bootstrap) = @_;
    chdir($self->{_tmpdir}); 
    my $output_base = $self->{_output_base};

    my $treeFile = $output_base."_fasttree.nwk";
    my $logFile = $output_base."_fasttree_log.txt";
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

    if ($bootstrap) {
        # use raxml to generate 100 data matrices
        # use fasttree to analyze them
        # use compareToBootstrap to count the per-branch support (as proportion)
        my $model = "GTRCAT";
        $model =  "PROTCATLG" if $self->{_alphabet} =~ /protein/i;
        print STDERR "alphabet = $self->{_alphabet} , model = $model\n" if $debug;
        @cmd = ('raxmlHPC-PTHREADS-SSE3',  '-f', 'j', '-b', '123', '-#', $bootstrap, '-s', $self->{_alignment_file}, '-n', 'boot_matrix', '-m', $model);
        $analysis_descriptor = "Generate bootstrap matrices";
        $self->add_analysis_step($analysis_descriptor, join(" ", @cmd));
        print STDOUT "run command: ". join(" ", @cmd)."\n" if $debug;
        my ($out, $err) = run_cmd(\@cmd);
        $self->add_analysis_out_err($out, $err);
        die "raxml failed to produce BS1" unless -f "$self->{_alignment_file}.BS1";
        my $multi_alignment_file = $self->{_alignment_file};
        $multi_alignment_file =~ s/\.{3,7}$//;
        $multi_alignment_file .= "_${bootstrap}BS.phy";
        @cmd = ("cat", "$self->{_alignment_file}.BS*", ">", $multi_alignment_file);
        print STDERR "Concatenate into one file:\n" . join(" ", @cmd) . "\n";
        system(join(" ", @cmd)); # requires shell interpolation
        system("rm $self->{_alignment_file}.BS*"); #allow shell expansion
        my $multi_tree_file = "multiple_bootstrap_trees.nwk";
        @cmd = ("FastTree", "-n", $bootstrap, "-out", $multi_tree_file, "-log", "fasttree_bootstrap_log");
        if ($self->{_alphabet} =~ /DNA/i) {
            push @cmd, "-nt", "-gtr";
        }
        elsif ($self->{_model} =~ /lg|wag/i) {
            push @cmd, "-" . lc($self->{_model});
        } # else defaults to jtt for proteins
        push @cmd, $multi_alignment_file;
        $analysis_descriptor = "Generate trees for each bootstrapped data matrix";
        $self->add_analysis_step($analysis_descriptor, join(" ", @cmd));
        print STDERR $analysis_descriptor, "\n", join(" ", @cmd), "\n" if $debug;
        print STDOUT "run command: ". join(" ", @cmd)."\n" if $debug;
        ($out, $err) = run_cmd(\@cmd);
        $self->add_analysis_out_err($out, $err);
        open CAT, ">>$logFile";
        open IN, "fasttree_bootstrap_log";
        print CAT <IN>;
        close CAT;
        $self->add_analysis_log($logFile);
        die "FastTree did not generate multiple tree file from boostrapped matrices" unless -f $multi_tree_file;
        my $tree_with_support = $output_base."_fasttree_boostrap_prop.nwk"; 
        @cmd = ("CompareToBootstrap", $treeFile, $multi_tree_file);
        $analysis_descriptor = "Map support onto best tree";
        $self->add_analysis_step($analysis_descriptor, join(" ", @cmd));
        print STDOUT "run command: ". join(" ", @cmd)."\n" if $debug;
        print STDOUT "run command: ". join(" ", @cmd)."\n" if $debug;
        ($out, $err) = run_cmd(\@cmd);
        $self->add_analysis_out_err($out, $err);
        open F, ">$tree_with_support";
        print F $out;
        close F;
        die "CompareToBootstrap did not generate tree with support values" unless -s $tree_with_support;
        $treeFile = $tree_with_support;
        return $treeFile;
    }

    move($treeFile, "$self->{_output_dir}/$treeFile");
    $self->add_analysis_tree($treeFile);
    $self->add_analysis_log($logFile);
    move($logFile, "$self->{_output_dir}/$logFile");
    chdir($self->{_original_wd}); 
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
