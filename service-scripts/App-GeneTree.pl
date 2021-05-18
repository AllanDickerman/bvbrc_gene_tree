#
# The GeneTree application.
use strict;
use Carp;
use Data::Dumper;
use File::Temp;
use File::Slurp;
use File::Basename;
use IPC::Run 'run';
use JSON;
use File::Copy ('copy', 'move');
use P3DataAPI;
use Bio::KBase::AppService::AppConfig;
use Bio::KBase::AppService::AppScript;
use Cwd;
use URI::Escape;
use Sequence_Alignment; # should be in lib directory
use Phylo_Tree; # should be in lib directory

our $global_ws;
our $global_token;

our $shock_cutoff = 10_000;

my $debug = 0;
$debug = $ENV{"GeneTreeDebug"} if exists $ENV{"GeneTreeDebug"};
Phylo_Tree::set_debug($debug); 
Phylo_Node::set_debug($debug); 
print STDERR "args = ", join("\n", @ARGV), "\n" if $debug;
our @analysis_step => ();# collect info on sequence of analysis steps
our @step_stack => (); # for nesting of child steps within parent steps

my $data_url = Bio::KBase::AppService::AppConfig->data_api_url;
#my $data_url = "http://www.alpha.patricbrc.org/api";

my $script = Bio::KBase::AppService::AppScript->new(\&build_tree, \&preflight);
my $rc = $script->run(\@ARGV);

sub preflight
{
    my($app, $app_def, $raw_params, $params) = @_;
    print STDERR "preflight: num params=", scalar keys %$params, "\n";
    my $pf = {
	cpu => 8,
	memory => "128G",
	runtime => 0,
	storage => 0,
	is_control_task => 0,
    };
    return $pf;
}

sub start_step {
    my $name = shift;
    my %step_details = {};
    push @analysis_step, \%step_details;
    $step_details{name} = $name;
    $step_details{comments} = ();
    $step_details{start_time} = time();
    print STDERR "start_step($step_details{name})\n";
    my $stack_depth = scalar @step_stack;
    my $current_step = $step_stack[$stack_depth-1];
    $step_details{parent_step} = $current_step;
    push @step_stack, $name;
    return \@{$step_details{comments}}, \%step_details
}

sub end_step {
    my $name = shift;
    print STDERR "end_step($name)\n";
    my $current_step = pop @step_stack;
    unless ($name eq $current_step) {
        print STDERR "Problem! at end_step name is wrong: $name should be $current_step";
    }
    my $step_index = scalar @analysis_step - 1;
    my $step_details = $analysis_step[$step_index];
    $step_details->{end_time} = time();
}

sub write_report {
    my ($output_file, $tree_graphic_file) = @_;
    print STDERR "write_report()\n";
    open F, ">$output_file";
    print F "<HTML>\n<h1>Gene Tree Output</h1>\n";
    if (-e $tree_graphic_file) {
        my $svg_text = read_file($tree_graphic_file);
        print F $svg_text, "\n\n";
    }
    print F "<h2>Analysis Steps</h2>\n";
    for my $step (@analysis_step) {
        print F "<h3>$step->{name}</h3>\n";
        my $time_string = localtime($step->{start_time});
        print F "<p>Start time $time_string\n";
        my $duration = $step->{end_time} - $step->{start_time};
        if ($duration > 10) { 
            print F "<p>Duration ", $duration, " seconds.\n";
        }
        if (scalar @{$step->{comments}}) {
            #print F "<p>Comments: ", scalar @{$step->{comments}}, "\n<ul>\n";
            print F "<ul>\n";
            for my $comment (@{$step->{comments}}) {
                print F "<li>$comment\n";
            }
            print F "</ul>\n";
        }
        if (exists $step->{stdout} and $step->{stdout} =~ /\S/) {
            print F "<div id=\"$step->{name}_stdout\"><pre>\n";
            print F $step->{stdout}, "\n</pre></div>\n";
        }
        if (exists $step->{stderr} and $step->{stderr} =~ /\S/) {
            print F "<div id=\"$step->{name}_stderr\"><pre>\n";
            print F $step->{stderr}, "\n</pre></div>\n";
        }
    }
    print F "</HTML>\n";
}

sub retrieve_sequence_data {
    my ($app, $params, $api, $tmpdir) = @_;
    my ($step_comments, $step_details) = start_step("retrieve_sequence_data");
    my $comment;
    my ($aligned_state, $any_in_memory) = (0, 0);
    $aligned_state = (scalar(@{$params->{sequences}}) == 1 and $params->{sequences}->[0]->{type} =~ /Aligned/i); 
    print STDERR "Number of sequence data sets = ", scalar(@{$params->{sequences}}), "\n";
    #push @{$step_comments}, "number of input sequence sets = " . scalar @{$params->{sequences}};
    for my $sequence_item (@{$params->{sequences}}) {
        print STDERR "data item: $sequence_item->{type}, $sequence_item->{filename}\n";
        push @{$step_comments}, "reading $sequence_item->{type} $sequence_item->{filename}";
        my $num_seqs = 0;
        if ($sequence_item->{type} =~ /FASTA/i) {
            # then it is one of the fasta formats in the workspace
            my $local_file = $sequence_item->{filename};
            $local_file =~ s/.*\///;
            print STDERR "About to dowload file $sequence_item->{filename} to $tmpdir/$local_file\n";
            $app->workspace->download_file($sequence_item->{filename}, "$tmpdir/$local_file", 1, $global_token);
            $sequence_item->{local_file} = $local_file;
            open F, "$tmpdir/$local_file";
            while (<F>) {
                $num_seqs++ if /^>/ }
            close F
        }
        elsif ($sequence_item->{type} eq "feature_group" or $sequence_item->{type} eq "feature_ids") {
            # need to get feature sequences from database 
            $any_in_memory = 1;
            $comment = "retrieving sequences for feature group $sequence_item->{filename}";
            push @{$step_comments}, $comment;
            my $feature_ids;
            if ($sequence_item->{type} eq 'feature_group') {
                print STDERR "retrieving ids for feature_group: $sequence_item->{filename}\n" if $debug;
                $feature_ids = $api->retrieve_patricids_from_feature_group($sequence_item->{filename});
            }
            else {
                $feature_ids = $sequence_item->{sequences}
            }
            print STDERR "feature_ids = ", join(", ", @$feature_ids), "\n" if $debug;
            if ($params->{alphabet} eq 'DNA') {
                $sequence_item->{sequences} = $api->retrieve_nucleotide_feature_sequence($feature_ids);
            }
            else {
                $sequence_item->{sequences} = $api->retrieve_protein_feature_sequence($feature_ids);
            }
            $num_seqs = scalar keys %{$sequence_item->{sequences}};
        }
        $comment = "Number of sequences retrieved = $num_seqs";
        push @{$step_comments}, $comment;
        print STDERR "$comment\n";
    }
    $comment = $aligned_state ? "sequences are aligned" : "sequences are unaligned";
    push @{$step_comments}, $comment;
    print STDERR "$comment\n";
    end_step("retrieve_sequence_data");
    return ($aligned_state, $any_in_memory);
}

sub merge_sequence_sets {
    my ($params, $tmpdir) = @_;
    my ($step_comments, $step_details) = start_step("merge_sequence_sets");
    my $comment = "merge_sequence_sets: $params, $tmpdir";
    push @{$step_comments}, $comment;
    print STDERR "$comment\n";
    # goal is to end up with all sequences written to a single unaligned FASTA file
    # only needed when multiple input sequence sets are passed
    for my $seq_item (@{$params->{sequences}}) {
        if (exists $seq_item->{local_file}) {
            $comment = "Reads sequences from $seq_item->{local_file}";
            push @{$step_comments}, $comment;
            print STDERR "$comment\n";
            open INDATA, "$tmpdir/$seq_item->{local_file}";
            my $seqid = '';
            $seq_item->{sequences} = \{};
            while (<INDATA>) {
                if (/^>(\S+)/) {
                    $seqid = $1;
                    if (exists $seq_item->{sequences}->{$seqid}) { # fix up non-unique seq identifier
                        $comment = "got non-unique sequence ID $seqid";
                        my $suffix = 2;
                        while (exists $seq_item->{sequences}->{$seqid . "_$suffix"}) {
                            $suffix++
                        }
                        $seqid = $seqid . "_$suffix";
                        $comment .= ", appending $suffix";
                        push @{$step_comments}, $comment;
                    }
                }
                else {
                    chomp;
                    $seq_item->{sequences}->{$seqid} .= $_
                }
            }
            $comment = "Found " . scalar keys %{$seq_item->{sequences}}, " sequences in $seq_item->{local_file}.";
            push @{$step_comments}, $comment;
            print STDERR "$comment\n";
        }
    }
    $comment = "merging all sequences";
    push @{$step_comments}, $comment;
    print STDERR "$comment\n";
    my %all_seqs = ();
    for my $seq_item (@{$params->{sequences}}) {
        $comment = "visiting $seq_item->{filename}";
        push @{$step_comments}, $comment;
        print STDERR "$comment\n";
        # each item should now have a hash of sequences
        for my $seqid (keys %{$seq_item->{sequences}}) {
            my $orig_seqid = $seqid;
            if (exists $all_seqs{$seqid}) {
                $comment = "got non-unique sequence ID $seqid";
                my $suffix = 2;
                while (exists $all_seqs{$seqid . "_$suffix"}) {
                    $suffix++
                }
                $seqid = $seqid . "_$suffix";
                $comment = ", appending $suffix";
                push @{$step_comments}, $comment;
            }
            $all_seqs{$seqid} = $seq_item->{sequences}->{$orig_seqid}
        }
    }
    $comment = "Total number of sequences is ". scalar keys %all_seqs;
    push @{$step_comments}, $comment;
    print STDERR "$comment\n";
    my $outfile = "$tmpdir/all_sequences_unaligned.fa";
    $comment = "write all sequences to $outfile";
    push @{$step_comments}, $comment;
    print STDERR "$comment\n";
    open F, ">$outfile";
    for my $seq_id (sort keys %all_seqs) {
        my $seq = $all_seqs{$seq_id};
        $seq =~ tr/-//d;
        $seq =~ tr/\s//d;
        print F ">$seq_id\n$seq\n";
    }
    close F;
    end_step("merge_sequence_sets");   
    return $outfile;
}

sub build_tree {
    my ($app, $app_def, $raw_params, $params) = @_;

    my $time1 = `date`;
    print "Proc GeneTree build_tree ", Dumper($app_def, $raw_params, $params);
    $global_token = $app->token()->token();
    # print STDERR "Global token = $global_token\n";
    my @outputs; # array of tuples of (filename, filetype)
    my $api = P3DataAPI->new;
    my $tmpdir = File::Temp->newdir( "/tmp/GeneTree_XXXXX", CLEANUP => !$debug );
    system("chmod", "755", "$tmpdir");
    print STDERR "created temp dir: $tmpdir\n";
   
    my ($aligned_state, $any_in_memory) = retrieve_sequence_data($app, $params, $api, $tmpdir);
    my $unaligned_file = undef;
    my $aligned_file = undef;
    if (scalar @{$params->{sequences}} > 1 or $any_in_memory) {
        #print STDERR "About to merge sequences, tmpdir=$tmpdir\n";
        $unaligned_file = merge_sequence_sets($params, $tmpdir)
    }
    elsif (!$aligned_state) {
        $unaligned_file = $params->{sequences}[0]->{local_file}
    }
    if ($unaligned_file) {
        $aligned_file = $unaligned_file;
        $aligned_file =~ s/\..{1,6}$//;
        $aligned_file =~ s/unaligned/aligned/;
        $aligned_file .= ".afa";
        run_muscle($unaligned_file, $aligned_file);
    }
    else { # if there was only one sequence source, and it was aligned and a file, leave it alone
        $aligned_file = "$tmpdir/$params->{sequences}[0]->{local_file}"
    }
    my $alignment = new Sequence_Alignment($aligned_file);
    my $seq_ids = $alignment->get_ids();
    print STDERR "seq ids: $seq_ids\n";
    print STDERR "seq ids: ", join(", ", @{$seq_ids}), "\n";
    #$params = localize_params($tmpdir, $params);
    #print "after localize_params:\n", Dumper($params);
    #
    #if ($debug) {
    #    push @outputs, [$jdesc, 'txt']
    #}

    #print STDERR "number of data inputs is " . scalar(@{$params->{sequences}}) . "\n";
    if ($params->{trim_threshold} or $params->{gap_threshold})
    {
        #$aligned_file = trim_alignment($aligned_file, $params, $tmpdir);
        my ($trim_comments, $trim_details) = start_step("trim_alignment");
        my $comment = "performing trimming on alignment";
        print STDERR "$comment\n";
        push @{$trim_comments}, $comment;
        $alignment->write_stats(*STDERR);
        if (exists $params->{trim_threshold} and $params->{trim_threshold} > 0) {
            $comment = "trim ends of alignment to density "+$params->{trim_threshold};
            print STDERR "$comment\n";
            push @{$trim_comments}, $comment;
            $alignment->end_trim($params->{trim_threshold});
            $alignment->write_stats(*STDERR);
        }
        if (exists $params->{gap_threshold} and $params->{gap_threshold} > 0) {
            $comment = "deleting sequences with gap proportion greater than "+$params->{gap_threshold};
            print STDERR "$comment\n";
            push @{$trim_comments}, $comment;
            $alignment->delete_gappy_seqs($params->{gap_threshold});
            $alignment->write_stats(*STDERR);
        }
        my $trimmed_alignment = $aligned_file;
        $trimmed_alignment =~ s/\.msa/_trimmed.msa/;
        $alignment->write_fasta("$tmpdir/$trimmed_alignment");

        my $data_type = "aligned_" . lc($params->{alphabet}) . "_fasta";
        push @outputs, ["$tmpdir/$trimmed_alignment", $data_type];
    }
    run("echo $tmpdir && ls -ltr $tmpdir");
    # count number of sequences
    my $num_seqs = $alignment->get_ntaxa();

    my $model = "AUTO"; # default for protein
    if ($params->{substitution_model}) {
        $model = $params->{substitution_model}
    }
    elsif ($params->{alphabet} =~ /DNA/i) {
        $model = "GTR"
    }
    my $output_name = $params->{output_file};
    my $alphabet = $params->{alphabet};
    my $recipe = "raxml"; #default
    if (defined $params->{recipe} and $params->{recipe}) {
        $recipe = lc($params->{recipe})
    }
    print STDERR "About to call tree program $recipe\n";
    if ($recipe eq 'raxml') {
        my $bad_chars = $alignment->write_fasta_for_raxml("$tmpdir/alignment_for_raxml.fa");
        my @tree_outputs = run_raxml("$tmpdir/alignment_for_raxml.fa", $alphabet, $model, $output_name, $tmpdir);
        if ($bad_chars) {
            for my $item (@outputs) {
                if ($item->[1] eq 'nwk') {
                    my $newick_file_name = "$tmpdir/$item->[0]";
                    open F, "+<", $newick_file_name;
                    my $raxml_newick = "";
                    while (<F>) {
                        chomp;
                        $raxml_newick .= $_
                    }
                    my $restored_newick = $alignment->restore_original_ids_in_raxml_tree($raxml_newick);
                    seek F, 0, 0;
                    truncate F, 0;
                    print F $restored_newick;
                    close F
                }

                # look for nkw file, revert any changes made to seq_ids for raxml
            }
        }

        push @outputs, @tree_outputs;
    } elsif ($recipe eq 'phyml') {
        $aligned_file =~ s/\.msa/.phy/;
        $alignment->write_phylip($aligned_file);
        my @tree_outputs = run_phyml($aligned_file, $alphabet, $model, $output_name, $tmpdir);
        push @outputs, @tree_outputs;
    } else {
        die "Unrecognized recipe: $recipe \n";
    }
    # generate tree graphic using figtree
    my $tree_file = undef;
    my $tree_graphic_file = undef;
    for my $output (@outputs) {
        my($ofile, $type) = @$output;
        if ($type eq 'nwk') {
            $tree_file = $ofile;
            $tree_graphic_file .= generate_tree_graphic($tree_file, $num_seqs);
        }
    }
    my $html_file = "$tmpdir/$params->{output_file}_gene_tree_report.html";
    write_report($html_file, $tree_graphic_file);
    push @outputs, [$html_file, "html"];
    
    if ($tree_file) {
        my $seq_ids_ar = $alignment->get_ids();
        print STDERR "ids: ", join(", ", @{$seq_ids_ar}), "\n";
        my $feature_metadata = get_feature_metadata($seq_ids_ar);
        my $tree = new Phylo_Tree($tree_file);
        my $phyloxml_data = $tree->write_phyloXML($feature_metadata);
        my $phyloxml_file = "$tmpdir/$params->{output_file}.xml";
        write_file($phyloxml_file, $phyloxml_data);
        push @outputs, [$phyloxml_file, "xml"];
    }
    
    print STDERR '\@outputs = '. Dumper(\@outputs);
    
    my $output_folder = $app->result_folder();
    # my $output_base   = $params->{output_file};
    #$app->workspace->create( { objects => [[$path, 'folder']] } );
    
    for my $output (@outputs) {
        my($ofile, $type) = @$output;
        next if $type eq 'folder';
        
        if (! -f $ofile) {
            warn "Output file '$ofile' of type '$type' does not exist\n";
            next;
        }
        
        my $filename = basename($ofile);
        #print STDERR "Output folder = $output_folder\n";
        print STDERR "Saving $filename => $output_folder as $type\n";
        $app->workspace->save_file_to_file($ofile, {}, "$output_folder/$filename", $type, 1,
                       (-s $ofile > $shock_cutoff ? 1 : 0), # use shock for larger files
                       $global_token);
    }
    my $time2 = `date`;
    write_output("Start: $time1"."End:   $time2", "$tmpdir/DONE");
}

sub trim_alignment {
    my ($aligned_file, $params, $tmpdir) = @_;
    my ($trim_comments, $details) = start_step("trim_alignment");
    my $comment = "performing trimming on alignment";
    print STDERR "$comment\n";
    push @{$trim_comments}, $comment;
    open my $ALIGNED, "$tmpdir/$aligned_file" or die "could not open aligned fasta";
    my $alignment = new Sequence_Alignment($ALIGNED);
    $alignment->write_stats(*STDERR);
    if (exists $params->{trim_threshold} and $params->{trim_threshold} > 0) {
        $comment = "trim ends of alignment to density "+$params->{trim_threshold};
        print STDERR "$comment\n";
        push @{$trim_comments}, $comment;
        $alignment->end_trim($params->{trim_threshold});
        $alignment->write_stats(*STDERR);
    }
    if (exists $params->{gap_threshold} and $params->{gap_threshold} > 0) {
        $comment = "deleting sequences with gap proportion greater than "+$params->{gap_threshold};
        print STDERR "$comment\n";
        push @{$trim_comments}, $comment;
        $alignment->delete_gappy_seqs($params->{gap_threshold});
        $alignment->write_stats(*STDERR);
    }
    my $trimmed_alignment = $aligned_file;
    $trimmed_alignment =~ s/\.msa/_trimmed.msa/;
    $alignment->write_fasta("$tmpdir/$trimmed_alignment");
    end_step("trim_alignment");
    return $trimmed_alignment;
}

sub run_raxml {
    my ($alignment_file, $alphabet, $model, $output_name, $tmpdir) = @_;
    my ($step_comments, $step_details) = start_step("run_raxml");
    print STDERR "In run_raxml, alignment = $alignment_file\n";
    my $parallel = $ENV{P3_ALLOCATED_CPU};
    $parallel = 2 if $parallel < 2;
    
    my $cwd = getcwd();
    
    if ($alphabet eq 'DNA') {
        $model = 'GTRGAMMA'
    }
    else {
        $model = "PROTCAT". uc($model)
    }

    my @cmd = ("raxmlHPC-PTHREADS-SSE3");
    push @cmd, ("-T", $parallel);
    push @cmd, ("-p", "12345");
    push @cmd, ("-m", $model);
    push @cmd, ("-s", $alignment_file);
    push @cmd, ("-n", $output_name);
    
    my $comment = "command = ". join(" ", @cmd);
    push @{$step_comments}, $comment;
    print STDERR "$comment\n\n";
   
    chdir($tmpdir); 
    my ($rc, $out, $err) = run_cmd(\@cmd);
    print STDERR "STDOUT:\n$out\n";
    print STDERR "STDERR:\n$err\n";
    $step_details->{stdout} = $out;
    $step_details->{stderr} = $err;
    
    my @outputs;
    my $bestTreeFile = $output_name . "_RAxML_bestTree.nwk";
    move("RAxML_bestTree.".$output_name, $bestTreeFile);
    push @outputs, ["$tmpdir/$bestTreeFile", 'nwk'];
    push @outputs, ["$tmpdir/RAxML_info.".$output_name, 'txt'];

    chdir($cwd);
    run("echo $tmpdir && ls -ltr $tmpdir");
    end_step("run_raxml");
    return @outputs;
}

sub run_phyml {
    my ($alignment_file, $alphabet, $model, $output_name, $tmpdir) = @_;
    my ($step_comments, $step_details) = start_step("run_phyml");

    my $cwd = getcwd();
    
    my $datatype = 'aa';
    if ($alphabet =~ /DNA/i) {
        $model = 'GTR';
        $datatype = 'nt'
    }

    my @cmd = ("phyml");
    push @cmd, ("-i", $alignment_file);
    push @cmd, ("-d", $datatype);
    push @cmd, ("-m", $model);
    
    print STDERR "cmd = ", join(" ", @cmd) . "\n\n";
   
    chdir($tmpdir); 
    my ($rc, $out, $err) = run_cmd(\@cmd);
    print STDERR "STDOUT:\n$out\n";
    print STDERR "STDERR:\n$err\n";
    $step_details->{stdout} = $out;
    $step_details->{stderr} = $err;
    
    my @outputs;
    my $treeFile = $alignment_file."_phyml_tree.nwk";
    move($alignment_file."_phyml_tree.txt", $treeFile);
    my $statsFile = $alignment_file."_phyml_stats.txt";
    push @outputs, [$treeFile, 'nwk'];
    push @outputs, [$statsFile, 'txt'];

    chdir($cwd);
    run("echo $tmpdir && ls -ltr $tmpdir");

    end_step("run_phyml");
    return @outputs;
}

sub get_feature_metadata {
    # return hashref from patric feature id to hash of field values
    my $ids = shift;
    print STDERR "in get_feature_metadata: ids: $ids\n";
    my ($step_comments, $step_details) = start_step("get_metadata");
    print STDERR "in get_feature_metadata: ids: ", join(", ", @{$ids}), "\n";
    print STDERR "junk\njunk\n";
    my $escaped = join(",", map { uri_escape $_ } @{$ids});
    my $url = "$data_url/genome_feature/?in(patric_id,($escaped))&select(patric_id,genome_id,product)";
    print STDERR "query=$data_url\n";
    my $resp = curl_text($url);
    print STDERR "response=$resp\n";
    if ( $resp =~ /Error/) {
        die "Problem! query did not return properly.\n";
    }
    my %id_hash = map { $_ => 1 } @$ids;
    my %feature_metadata = ();
    my %genome_metadata = ();
    $resp =~ tr/"//d;
    for my $member (split("\n", $resp)) {
        my ($id, $genome_id, $product) = split("\t", $member);
        if (exists $id_hash{$id}) {
            $feature_metadata{$id} = {};
            $feature_metadata{$id}{'genome_id'} = $genome_id;
            $feature_metadata{$id}{'product'} = $product;
            $genome_metadata{$genome_id} = {};
            print STDERR "added genome_id $genome_id for $id\n";
        }
    }
    #print STDERR "genome keys = ", join(", ", keys %genome_metadata), "\n\n";
    $escaped = join(",", keys %genome_metadata);
    $url = "$data_url/genome/?in(genome_id,($escaped))&select(genome_id,genus,species,host_name,geographic_location,collection_date)";
    print STDERR "query=$data_url\n";
    $resp = curl_text($url);
    print STDERR "response=$resp\n";
    $resp =~ tr/"//d;
    for my $member (split "\n", $resp) {
        my ($genome_id, $genus, $species, $host, $geography, $date) = split("\t", $member);
        %{$genome_metadata{$genome_id}} = ("taxon" => $genus."_".$species, "host" => $host, "geography" => $geography, 'date' =>$date);
        print STDERR "try this: $genome_id, $genome_metadata{$genome_id}\n"; 
        print STDERR "keys of gm{$genome_id}: ", join(", ", keys %{$genome_metadata{$genome_id}}), "\n";
    }
    for my $feature_id (keys %feature_metadata) {
        my $genome_id = $feature_metadata{$feature_id}{'genome_id'};
        print STDERR "feature $feature_id, genome $genome_id:\n ";
        for my $key (keys %{$genome_metadata{$genome_id}}) {
            print STDERR "\tset $key ";
            print STDERR "to $genome_metadata{$genome_id}->{$key}\n";
            $feature_metadata{$feature_id}->{$key} = $genome_metadata{$genome_id}->{$key}
        }
    }
    end_step("get_metadata");
    return \%feature_metadata;
}

sub get_md5_for_feature_group {
    #" return hashref from md5 to array of patric feature id(s) - one-to-many
    my ($group, $seq_type) = @_;
    my $seq_field_md5 = (lc($seq_type) eq "protein") ? 'aa_sequence_md5' : 'na_sequence_md5';
    my $escaped = uri_escape($group);
    my $url = "$data_url/genome_feature/?in(feature_id,FeatureGroup($escaped))&select(patric_id,$seq_field_md5)&http_accept=application/json&limit(25000)";
    my $resp = curl_json($url);
    my %md5_to_patric_id;
    for my $member (@$resp) {
        my $md5 = $member->{$seq_field_md5};
        $md5_to_patric_id{$md5} = () unless exists($md5_to_patric_id{$md5});
        push @{$md5_to_patric_id{$md5}}, $member->{'patric_id'}
    }
    if ($debug)
    {
        print "Number of patric IDs from group: ", scalar %md5_to_patric_id, "\n";
        for my $i (0..5) {
            print "member $i: $resp->[$i]{'patric_id'}, $resp->[$i]{aa_sequence_md5}\n"; 
        }
    }
    return \%md5_to_patric_id;
}

sub get_md5_for_feature_ids {
    #" return hashref from md5 to patric feature id
    my ($ids, $seq_type) = @_;
    my $seq_field_md5 = (lc($seq_type) eq "protein") ? 'aa_sequence_md5' : 'na_sequence_md5';
    my $escaped = join(",", map { uri_escape $_ } @$ids);
    my $url = "$data_url/genome_feature/?in(patric_id,($escaped))&select(patric_id,$seq_field_md5)&http_accept=application/json&limit(25000)";
    my $resp = curl_json($url);
    my %md5_to_patric_id;
    for my $member (@$resp) {
        my $md5 = $member->{$seq_field_md5};
        $md5_to_patric_id{$md5} = () unless exists($md5_to_patric_id{$md5});
        push @{$md5_to_patric_id{$md5}}, $member->{'patric_id'}
    }
    if ($debug)
    {
        print "Number of patric IDs from group: ", scalar %md5_to_patric_id, "\n";
        for my $i (0..5) {
            print "member $i: $resp->[$i]{'patric_id'}, $resp->[$i]{aa_sequence_md5}\n"; 
        }
    }
    return \%md5_to_patric_id;
}

sub get_feature_sequences_by_md5 {
    my $md5_to_patric_id = shift;
    my $list = join(",", keys(%$md5_to_patric_id));
    my $url = "$data_url/feature_sequence/?in(md5,($list))&select(md5,sequence)&http_accept=application/json&limit(25000)";
    my $resp = curl_json($url);
    print "response = \n\n$resp\n";
    my %patric_id_to_sequence;
    my $i = 0;
    for my $member (@$resp) {
        my $md5 = $member->{md5};
        my $sequence = $member->{'sequence'};
        for my $patric_id (@{$md5_to_patric_id->{$md5}})
        { # handle potential one-to-many relationship 
            $patric_id_to_sequence{$patric_id} = $sequence;
        }
        if ($debug and $i < 5) {
            print STDERR "$md5\t", join(",", @{$md5_to_patric_id->{$md5}}), "\t$sequence\n";
            $i++;
        }
    }
    return \%patric_id_to_sequence;
}

sub generate_tree_graphic {
    my ($input_newick, $num_tips) = @_;
    my ($step_comments, $step_details) = start_step("generate_tree_graphic");
    my $file_base = $input_newick;
    $file_base =~ s/\..{2,6}//;
    my $tree_graphic_file = "$file_base.svg";
    my $nexus_file = "$file_base.nex";
    my $comment = "run figtree input = $input_newick, output = $tree_graphic_file";
    push @{$step_comments}, $comment;
    print STDERR "$comment\n";
    open F, ">$nexus_file";
    print F "#NEXUS\nbegin trees;\n";
    my $tree_text = read_file($input_newick);
    print F "tree one = [&U] $tree_text\nend;\n\n";
    print F "begin figtree;\n";
    print F "set appearance.branchLineWidth=3.0;\n";
    print F "set tipLabels.fontSize=14;\n";
    print F "set tipLabels.fontName=\"sansserif\";\n";
    print F 'set trees.order=true;';
    print F 'set trees.orderType="increasing";';
    print F 'set trees.rooting=true;';
    print F 'set trees.rootingType="Midpoint";';
    print F "end;\n";
    close F;

    my @cmd = ("figtree", "-graphic", 'SVG');
    
    if ($num_tips > 40) {
        my $height = 600 + 15 * ($num_tips - 40); # this is an empirical correction factor to avoid taxon name overlap
        push @cmd, '-height', $height;
    }
    push @cmd, $nexus_file, $tree_graphic_file;
    my ($rc, $stdout, $stderr) =  run_cmd(\@cmd);
    $step_details->{stdout} = $stdout;
    $step_details->{stderr} = $stderr;
    end_step("generate_tree_graphic");
    return $tree_graphic_file
}

sub run_muscle {
    my ($unaligned, $aligned) = @_;
    print STDERR "run_muscle($unaligned, $aligned)\n";
    my ($step_comments, $step_details) = start_step("run_muscle");
    my $comment = "input = $unaligned, output = $aligned";
    push @{$step_comments}, $comment;
    print STDERR "$comment\n";
    my $cmd = ["muscle", "-in", $unaligned, "-out", $aligned];
    my ($rc, $stdout, $stderr) =  run_cmd($cmd);
    end_step("run_muscle");
    $step_details->{stdout} = $stdout;
    $step_details->{stderr} = $stderr;
}

sub curl_text {
    my ($url) = @_;
    my @cmd = ("curl", $url, curl_options());
    my $cmd = join(" ", @cmd);
    #$cmd =~ s/sig=[a-z0-9]*/sig=XXXX/g;
    print STDERR "$cmd\n";
    my ($out) = run_cmd(\@cmd);
    return $out;
}

sub curl_json {
    my ($url) = @_;
    my $out = curl_text($url);
    my $hash = JSON::decode_json($out);
    return $hash;
}

sub curl_options {
    my @opts;
    my $token = $global_token;
    push(@opts, "-H", "Authorization:$token");
    push(@opts, "-H", "Accept:text/tsv");
    #push(@opts, "-H", "Content-Type: multipart/form-data");
    return @opts;
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

sub get_ws {
    return $global_ws;
}

sub get_token {
    return $global_token;
}

sub get_ws_file {
    my ($tmpdir, $id) = @_;
    # return $id; # testing
    my $ws = get_ws();
    my $token = get_token();
    
    my $base = basename($id);
    my $file = "$tmpdir/$base";
    # return $file; # testing
    
    my $fh;
    open($fh, ">", $file) or die "Cannot open $file for writing: $!";

    print STDERR "GET WS => $tmpdir $base $id\n";
    system("ls -la $tmpdir");

    eval {
	$ws->copy_files_to_handles(1, $token, [[$id, $fh]]);
    };
    if ($@)
    {
	die "ERROR getting file $id\n$@\n";
    }
    close($fh);
    print "$id $file:\n";
    system("ls -la $tmpdir");

    return $file;
}

sub write_output {
    my ($string, $ofile) = @_;
    open(F, ">$ofile") or die "Could not open $ofile";
    print F $string;
    close(F);
}

sub verify_cmd {
    my ($cmd) = @_;
    system("which $cmd >/dev/null") == 0 or die "Command not found: $cmd\n";
}
