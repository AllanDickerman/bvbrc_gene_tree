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
use List::Util ('any');
use P3DataAPI;
use Bio::KBase::AppService::AppConfig;
use Bio::KBase::AppService::AppScript;
use Cwd;
use URI::Escape;
use Sequence_Alignment; # should be in lib directory
use Tree_Builder;

our $global_ws;
our $global_token;

our $shock_cutoff = 10_000;
my $max_genome_length = 250_000; #most/all single-sequence viruses are less than this
my @default_genome_metadata_fields = (
        "genome_name", "strain", "genbank_accessions", "subtype", "lineage", "clade", "h1_clade_global", "h1_clade_us", "h5_clade", "host_group", "host_common_name", "host_scientific_name", "collection_year", "geographic_group", "isolation_country", "state_province");
my @default_feature_metadata_fields = ("product", "accession", "patric_id");

our $debug = 0;
$debug = $ENV{"GeneTreeDebug"} if exists $ENV{"GeneTreeDebug"};
if ($debug) {
    print STDERR "debug = $debug\n" if $debug;
    Tree_Builder::set_debug($debug);
    Sequence_Alignment::set_debug($debug);
    print STDERR "args = ", join("\n", @ARGV), "\n";
}
our @analysis_step => ();# collect info on sequence of analysis steps
our @step_stack => (); # for nesting of child steps within parent steps
#my @original_sequence_ids; # list of items requested, can be different from those actually obtained
#my $sequence_identifier_type; # feature_id or genome_id or user_specified

my $data_url = Bio::KBase::AppService::AppConfig->data_api_url;
#$data_url = "https://patricbrc.org/api" if $debug;
print STDERR "data_url=\n$data_url\n" if $debug;

my $script = Bio::KBase::AppService::AppScript->new(\&build_tree, \&preflight);
my $rc = $script->run(\@ARGV);

sub preflight
{
    my($app, $app_def, $raw_params, $params) = @_;
    print STDERR "preflight: num params=", scalar keys %$params, "\n";
    my $api = P3DataAPI->new;
    print STDERR "Number of sequence data sets = ", scalar(@{$params->{sequences}}), "\n";
    my $num_seqs = 0;
    my $total_length = 0;
    for my $sequence_source (@{$params->{sequences}}) {
        print STDERR "data item: $sequence_source->{type}, $sequence_source->{filename}\n";
        if ($sequence_source->{type} =~ /FASTA/i) {
            # then it is one of the fasta formats in the workspace
            my $fasta_string = $app->workspace->download_file_to_string($sequence_source->{filename}, $global_token);
            for $_ (split("\n", $fasta_string)) {
                if (/^>/) {
                    $num_seqs++;
                }
                else {
                    $total_length += length($_);
                }
            }
        }
        elsif ($sequence_source->{type} eq "feature_group") {
            # need to get feature sequences from database 
            my $feature_group = $sequence_source->{filename};
            my $na_or_aa = ('aa', 'na')[$params->{alphabet} eq 'DNA']; 
            my $seq_list = $api->retrieve_sequences_from_feature_group($feature_group, $na_or_aa);
            for my $item (@$seq_list) {
                $num_seqs++;
                $total_length += length($item->{sequence});
            }
        }
        elsif ($sequence_source->{type} eq "genome_group") {
            my $genome_group = $sequence_source->{filename};
            my $genome_ids = $api->retrieve_patric_ids_from_genome_group($genome_group);
            $num_seqs += scalar @{$genome_ids};
            #print STDERR "Num seqs = $num_seqs, genome_ids = @{$genome_ids}\n";
            my @genome_validation_data = $api->retrieve_genome_metadata($genome_ids, ['genome_length']);
            for my $info (@genome_validation_data) {
                $total_length += $info->{genome_length};
                print STDERR "$total_length\n";
            }
        }
        print STDERR "Num seqs = $num_seqs, total lenth = $total_length\n";
    }
    my $run_time = int($total_length/ 200);
    $run_time = 600 if $run_time < 600; # give a minimum of 10 minutes
    print STDERR "caclulating runtime as total_length/200 = $total_length / 200 = $run_time\n";
    my $pf = {
	cpu => 8,
	memory => "32G",
	runtime => $run_time,
	storage => 0,
	is_control_task => 0,
    };
    return $pf;
}

sub start_step {
    my $name = shift;
    my %step_info = {};
    push @analysis_step, \%step_info;
    $step_info{name} = $name;
    $step_info{comments} = ();
    $step_info{start_time} = time();
    print STDERR "start_step($step_info{name})\n";
    my $stack_depth = scalar @step_stack;
    my $current_step = $step_stack[$stack_depth-1];
    $step_info{parent_step} = $current_step;
    push @step_stack, $name;
    return \@{$step_info{comments}}, \%step_info
}
sub add_step_command_line {
    my $cmd = shift;
    my $last_index = $#analysis_step;
    my $step_info = $analysis_step[$last_index];
    $step_info->{command_line} = $cmd;
}

sub end_step {
    my $name = shift;
    print STDERR "end_step($name)\n";
    my $current_step = pop @step_stack;
    #unless ($name eq $current_step) {
    #    print STDERR "Problem! at end_step name is wrong: $name should be $current_step" if $debug;
    #}
    my $step_index = scalar @analysis_step - 1;
    my $step_info = $analysis_step[$step_index];
    $step_info->{end_time} = time();
}

sub add_analysis_step { #allow adding a step recorded by Tree_Builder object
    my $step_info = shift;
    push @analysis_step, $step_info;
}

sub build_sequence_length_table {
    # $lengths is a hashref: {id}{locus} = length, as returned by Sequence_Alignment::get_sequence_lengths
    # $loci is an arrayref of all locus names known for the alignment (used to decide whether to show per-locus columns)
    # one row per Sequence ID; one column per locus, unless there is only one locus (then a single Length column)
    my ($lengths, $loci) = @_;
    my $multi_locus = (scalar @$loci > 1);
    my $html = "<table border=\"1\" cellpadding=\"3\" cellspacing=\"0\">\n<tr><th>Sequence ID</th>";
    if ($multi_locus) {
        $html .= "<th>$_</th>" for @$loci;
    }
    else {
        $html .= "<th>Length</th>";
    }
    $html .= "</tr>\n";
    for my $id (sort keys %$lengths) {
        $html .= "<tr><td>$id</td>";
        if ($multi_locus) {
            for my $locus (@$loci) {
                my $length = $lengths->{$id}{$locus};
                $html .= "<td>" . (defined $length ? $length : '-') . "</td>";
            }
        }
        else {
            my ($length) = values %{$lengths->{$id}};
            $html .= "<td>" . (defined $length ? $length : '-') . "</td>";
        }
        $html .= "</tr>\n";
    }
    $html .= "</table>\n";
    return $html;
}

sub add_sequence_summary_step {
    # record a report step consisting only of a sequence-length table (no command line)
    my ($name, $lengths, $loci) = @_;
    my $now = time();
    add_analysis_step({
        name => $name,
        start_time => $now,
        end_time => $now,
        table => build_sequence_length_table($lengths, $loci),
    });
}

sub build_alignment_occupancy_table {
    # $occupancy is a hashref: {locus} => { length => alignment length, percent_non_gap => ... }
    # as returned by Sequence_Alignment::get_alignment_occupancy
    # $loci is an arrayref of locus names, in the order rows should be listed
    # one row per locus; Locus column is omitted when there is only one locus
    my ($occupancy, $loci) = @_;
    my $multi_locus = (scalar @$loci > 1);
    my $html = "<table border=\"1\" cellpadding=\"3\" cellspacing=\"0\">\n<tr>";
    $html .= "<th>Locus</th>" if $multi_locus;
    $html .= "<th>Alignment Length</th><th>% Non-gap</th></tr>\n";
    for my $locus (@$loci) {
        my $length = $occupancy->{$locus}{length};
        my $percent = $occupancy->{$locus}{percent_non_gap};
        $html .= "<tr>";
        $html .= "<td>$locus</td>" if $multi_locus;
        $html .= "<td>$length</td><td>" . sprintf("%.1f", $percent) . "%</td></tr>\n";
    }
    $html .= "</table>\n";
    return $html;
}

sub add_alignment_occupancy_step {
    # record a report step consisting only of a per-locus alignment occupancy table (no command line)
    my ($name, $occupancy, $loci) = @_;
    my $now = time();
    add_analysis_step({
        name => $name,
        start_time => $now,
        end_time => $now,
        table => build_alignment_occupancy_table($occupancy, $loci),
    });
}

sub write_report {
    my ($output_file, $title, $tree_graphic_files) = @_;
    print STDERR "write_report()\n";
    open F, ">$output_file";
    print F "<HTML>\n<head>\n";
    print F "<script>\nfunction toggle_visibility(element_name) {
      var x = document.getElementById(element_name);
        if (x.style.display == \"none\") {
              x.style.display = \"block\";
          } else {
              x.style.display = \"none\";
          }
      }\n</script>\n";
    print F "</head><body>\n<h1>$title</h1>\n";
    for my $tree_graphic_file (@$tree_graphic_files) {
        if (-e $tree_graphic_file) {
            my $element_name = 'tree_plot_' . $tree_graphic_file;
            print F "FigTree Plot $tree_graphic_file: <button onclick=\"toggle_visibility('$element_name')\">Show/Hide</button>\n";
            print F "<div id=\"$element_name\" style=\"display:block; background:#ffffff\" \n";
            print F "    onclick=\"toggle_visibility('$element_name')\">\n";
            my $svg_text = read_file($tree_graphic_file);
            print F $svg_text, "\n</div><br>\n";
        }
        else { print STDERR "Tree graphic file not found: $tree_graphic_file\n"; }
    }
    print F "<h2>Analysis Steps</h2>\n";
    for my $step (@analysis_step) {
        print F "<h3>$step->{name}</h3>\n";
        if (exists $step->{table}) {
            print F $step->{table}, "\n";
        }
        if (exists $step->{command_line}) {
            print F "<pre>$step->{command_line}</pre>\n";
        }
        if (exists $step->{comments} and scalar @{$step->{comments}}) {
            print F "Annotation fields:\n" if $step->{name} =~ /Write PhyloXML/;
            print F "<ul>\n";
            for my $comment (@{$step->{comments}}) {
                print F "<li>$comment\n";
            }
            print F "</ul>\n";
        }
        if (exists $step->{details} and $step->{details} =~ /\S/) {
            my $element_name = "$step->{name}_details";
            $element_name =~ tr/ /_/;
            print F "Details: <button onclick=\"toggle_visibility('$element_name')\">Show/Hide</button>\n";
            print F "<div id=\"$element_name\" style=\"display:none; background:#f0f0f0\" \n";
            print F "    onclick=\"toggle_visibility('$element_name')\">\n";
            print F "<pre>\n", $step->{details}, "\n</pre></div>\n";
        }
        my $duration = $step->{end_time} - $step->{start_time};
        if ($duration > 10) { 
            print F "<p>Duration ", $duration, " seconds.\n";
        }
    }
    my $start_time = $analysis_step[0]->{start_time};
    my $time_string = localtime($start_time);
    my $duration = time() - $start_time;
    print F "<p>Start time $time_string<br>\n";
    print F "Duration: $duration<br>\n";
    print F "</body></HTML>\n";
}

sub validate_genomes_for_viral_tree {
    my($genome_ids, $api) = @_;
    print STDERR "examine genome metadata to test for single-sequence virus under $max_genome_length:\n";
    print STDERR "genome_ids = ", join(", ", @$genome_ids), ".\n";
    my $genome_validation_fields = ['genome_id', 'contigs', 'superkingdom', 'genome_length'];
    my @genome_validation_data = $api->retrieve_genome_metadata($genome_ids, $genome_validation_fields);
    print STDERR join("\t", @$genome_validation_fields), "\n";
    for my $info (@genome_validation_data) {
        for my $key ('genome_id', 'contigs', 'superkingdom', 'genome_length') {
            print STDERR "$info->{$key}\t";
        }
        print STDERR "\n";
        if ($info->{'genome_length'} > $max_genome_length) {
            print STDERR "Problem: length of genome $info->{'genome_id'} exceeds $max_genome_length ($info->{'genome_length'}).\nExiting.\n";
            `cd ..`;
            exit(1);
        }
        if ($info->{'superkingdom'} ne 'Viruses') {
            print STDERR "Problem: superkingdom of genome $info->{'genome_id'} is not 'Viruses' ($info->{'superkingdom'}).\nExiting.\n";
            `cd ..`;
            exit(1);
        }
        if ($info->{'contigs'} > 1) {
            print STDERR "Genome $info->{'genome_id'} has multiple contigs (segments) ($info->{'contigs'}).\n";
        }
    }
    print STDERR "All genomes are viruses and all are under $max_genome_length bases.\n";
}


sub retrieve_sequence_data {
    my ($app, $params, $api) = @_;
    my ($step_comments, $step_info) = start_step("Gather Sequence Data");
    my $seq_al = new Sequence_Alignment(); # object to store sequences
    my $comment;
    my ($aligned_state, $any_in_memory) = (0, 0);
    $aligned_state = (scalar(@{$params->{sequences}}) == 1 and $params->{sequences}->[0]->{type} =~ /Aligned/i); 
    print STDERR "Number of sequence data sets = ", scalar(@{$params->{sequences}}), "\n";
    my $total_seqs = 0;
    my @empty_sequences; # keep track of entries lacking sequence data
    for my $sequence_source (@{$params->{sequences}}) {
        $comment = "fetch $sequence_source->{type} $sequence_source->{filename}";
        print STDERR $comment;
        push @{$step_comments}, $comment;
        my $num_seqs = 0;
        if ($sequence_source->{type} =~ /FASTA/i) {
            # then it is one of the fasta formats in the workspace
            my $local_file = $sequence_source->{filename};
            $local_file =~ s/.*\///; #remove path, leave filename
            print STDERR "About to copy file $sequence_source->{filename} to $local_file\n";
            $app->workspace->download_file($sequence_source->{filename}, $local_file, 1, $global_token);
            my @ids = $seq_al->read_file($local_file);
            $num_seqs = scalar @ids;
            for my $id (@ids) {
                $seq_al->set_metadata($id, 'source', $local_file); #mark where this data came from
                if ($id =~ /^fig\|\d+\.\d+\..{3}\.\d+$/) {
                    print "try user identifier as patric_id: $id\n" if $debug;
                    $seq_al->set_metadata($id, 'database_link', 'patric_id');
                }
                elsif ($id =~ /^\w+\.\d+\.\d+\.\w+\.\w+\.\d+\.\d+\.(fwd|rev)$/) {
                    print "try user identifier as feature_id: $id\n" if $debug;
                    $seq_al->set_metadata($id, 'database_link', 'feature_id');
                }
                elsif ($id =~ /(\d+\.\d+$)/) {
                    print "try user identifier as genome_id $id\n" if $debug;
                    $seq_al->set_metadata($id, 'database_link', 'genome_id');
                }
                elsif ($id =~ /(acc\|):?(\d+\.\d+$)/) {
                    print "try user identifier as genome_id $id\n" if $debug;
                    $seq_al->set_metadata($id, 'database_link', 'genome_id');
                }
            }
        }
        elsif ($sequence_source->{type} eq "feature_group") {
            # need to get feature sequences from database 
            my $feature_group = $sequence_source->{filename};
            if ($debug) {
                my $feature_ids = $api->retrieve_patricids_from_feature_group($feature_group);
                print STDERR "\nfeature_ids = @$feature_ids\n";
                print STDERR "number feature_ids is ", scalar @$feature_ids, "\n";
            }

            my $seq_list;
            my $na_or_aa = ('aa', 'na')[$params->{alphabet} eq 'DNA']; 
            $seq_list = $api->retrieve_sequences_from_feature_group($feature_group, $na_or_aa);
            if ($debug) {
                print STDERR "retreive_sequences_from_feature_group return is $seq_list\n"; 
                print STDERR "length is ", scalar @$seq_list, "\n"; 
                print STDERR "First element = $seq_list->[0]\n";
                for my $key (sort keys %{$seq_list->[0]}) {
                    print STDERR "    \t$key\t$seq_list->[0]->{$key}\n";
                }
            }
            $feature_group =~ s/.*\///; # remove path preceding name of feature group
            for my $item (@$seq_list) {
               my $id = $seq_al->add_seq($item->{feature_id}, $item->{sequence});
               $seq_al->set_metadata($id, 'data_source', $feature_group);
               $seq_al->set_metadata($id, 'database_link', 'feature_id');
            }
            $num_seqs = scalar @$seq_list;
            $comment = "number of sequence fetures retrieved from $feature_group: $num_seqs";
            print STDERR $comment, "\n";
        }
        elsif ($sequence_source->{type} eq "genome_group") {
            my $genome_group = $sequence_source->{filename};
            my $genome_ids = $api->retrieve_patric_ids_from_genome_group($genome_group);
            print STDERR "got genome ids: "+join(", ", @$genome_ids), ".\n";
            #  do we need to validate length and that each is a virus? 
            #validate_genomes_for_viral_tree($genome_ids, $api);

            my @segments_to_use;
            if ($params->{genome_selection}{selected_segments}) {
                print STDERR "Limit analysis to segments: ", join(", ", @{$params->{genome_selection}{selected_segments}}), "\n";
                @segments_to_use = @{$params->{genome_selection}{selected_segments}};
            }
            $genome_group =~ s/.*\///; # remove path preceding name of genome group
            for my $genome_id (@$genome_ids) {
                my ($resp, $data) = $api->submit_query('genome_sequence', "eq(genome_id,$genome_id)", "sequence,contig");
                #print "for $genome_id: resp = $resp\tdata=$data\tdata->[0]=$data->[0]\n" if $debug;
                for my $record (@$data) {
                    my $use = 1;
                    my $locus = undef;
                    if ($record->{segment}) {
                        if (@segments_to_use) {
                            $use = any { $_ eq $record->{segment} } @segments_to_use;
                        }
                        $locus = "segment_" . $record->{segment};
                    }
                    if ($use) {
                        my $id = $seq_al->add_seq($genome_id, $record->{sequence}, $locus);
                        $seq_al->set_metadata($id, "data_source", $genome_group);
                        $seq_al->set_metadata($id, "database_link", "genome_id");
                    }
                    else {
                        print STDERR "skipping segment $record->{segment} not in segments_to_use\n";
                    }
                }
            }
            $num_seqs = scalar(@$genome_ids);
        }
        elsif ($sequence_source->{type} eq "feature_ids") {
            # need to get feature sequences from database 
            my $feature_seq = $api->retrieve_protein_feature_sequence($sequence_source->{feature_ids});
            for my $patric_id (keys %$feature_seq) {
                my $id = $seq_al->add_seq($patric_id, $feature_seq->{$patric_id}, "feature_list");
                $seq_al->set_metadata($id, "database_link", "feature_id");
                $seq_al->set_metadata($id, "data_source", "feature_list");
            }
            $num_seqs = scalar keys %{$feature_seq};
        }
        $comment = "$num_seqs entries retrieved from $sequence_source->{type} $sequence_source->{filename}\n";
        push @{$step_comments}, $comment;
        print STDERR $comment;
    }
    if (scalar @{$params->{sequences}} > 1) {
        my $num_seqs = $seq_al->get_ntaxa();
        $comment = "total sequnces retrieved = $num_seqs\n";
        push @{$step_comments}, $comment;
        print STDERR "$comment\n";
    }
    if (scalar @empty_sequences) {
        $comment = "records lacking sequence data: " . join(", ", @empty_sequences);
        push @{$step_comments}, $comment;
        print STDERR "$comment\n";
    }
    end_step("Gather Sequence Data");
    return $seq_al;
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
    print STDERR "created temp dir: $tmpdir, cleanup = ", !$debug, "\n";
    my $original_wd = getcwd();
    chdir($tmpdir); # do all work in temporary 
   
    my @feature_metadata_fields = @default_feature_metadata_fields;
    if (exists $params->{feature_metadata_fields}) {
        @feature_metadata_fields = @{$params->{feature_metadata_fields}};
    }
    #ensure that feature_id and genome_id are retrieved
    push @feature_metadata_fields, "feature_id" unless any { $_ eq 'feature_id'} @feature_metadata_fields;
    push @feature_metadata_fields, "genome_id" unless any { $_ eq 'genome_id'} @feature_metadata_fields;

    my @genome_metadata_fields = @default_genome_metadata_fields;
    if (exists $params->{genome_metadata_fields}) {
        @genome_metadata_fields = @{$params->{genome_metadata_fields}};
    }

    #my $seq_list = retrieve_sequence_data($app, $params, $api);
    my $seq_al = retrieve_sequence_data($app, $params, $api);
    my $num_seqs = $seq_al->get_ntaxa();
    print STDERR "After retrieval, number of sequences is $num_seqs\n" if $debug;
    if ($num_seqs < 4) { #need at least 4 seuqences to build a tree
        print STDERR "After retrieval, number of sequences is $num_seqs, less than 4. Cannot build a tree.\n";
        exit(1);
    }
    my @loci = $seq_al->get_locus_ids();
    add_sequence_summary_step("Sequence Summary (Before Alignment)", $seq_al->get_sequence_lengths(), \@loci);

    my $is_aligned = 0;
    my $outfile;
    if (scalar @{$params->{sequences}} == 1 and $params->{sequences}->[0]->{type} =~ /Aligned/i) {
        $is_aligned = 1;
    }
    else {
        my ($step_comments, $step_info) = start_step("Align Sequences");
        my ($cmd_lines, $stdout) = $seq_al->align();
        $step_info->{"command_line"} = join("", @{$cmd_lines});
        $step_info->{"stdout"} = $stdout;
        end_step();
    }
    if ($params->{trim_threshold})
    {
        my ($step_comments, $step_info) = start_step("End-Trim Alignment");
        my ($cmd_lines, $stdout) = $seq_al->end_trim($params->{trim_threshold});
        $step_info->{"command_line"} = join("", @$cmd_lines);
        $step_info->{"stdout"} = $stdout;
        end_step();
    }
    if ($params->{gap_threshold})
    {
        my ($step_comments, $step_info) = start_step("Filter Gappy Seqs");
        my ($cmd_lines, $stdout) = $seq_al->delete_gappy_seqs($params->{gap_threshold});
        $step_info->{"command_line"} = join("", @$cmd_lines);
        $step_info->{"stdout"} = $stdout if $stdout;
        end_step();
    }
    add_alignment_occupancy_step("Alignment Occupancy", $seq_al->get_alignment_occupancy(), \@loci);

    unless ($params->{recipe}) {
        $params->{recipe} = 'fasttree';
    }
    my $alphabet = $params->{alphabet};
    run("echo $tmpdir && ls -ltr $tmpdir") if $debug;

    my $model = "LG"; # default for protein
    if ($params->{substitution_model}) {
        $model = $params->{substitution_model}
    }
    elsif ($params->{alphabet} =~ /DNA/i) {
        $model = "GTR"
    }
    my @tree_outputs;

    my $threads = 2;
    if (exists $ENV{P3_ALLOCATED_CPU}) {
        $threads = $ENV{P3_ALLOCATED_CPU};
        print STDERR "P3_ALLOCATED_CPU = $ENV{P3_ALLOCATED_CPU}\n";
    }
    print STDERR "Tree program = $params->{recipe}\n";

    my @segments_to_tree = $seq_al->get_locus_ids();
    # add the empty locus if we intend to concatenate
    if ($params->{genome_selection} and $params->{genome_selection}{concat_segments}) {
        push(@segments_to_tree, "concatenated");
    }
    print STDERR "segments to build trees: " . join(" ", @segments_to_tree) . "\n" if $debug;

    my ($step_comments, $step_info) = start_step("Build Tree using $params->{recipe}");
    for my $alignment_component (@segments_to_tree) {
        print STDERR "build tree for $alignment_component\n" if $debug;
        my $alignment_file_base = $params->{output_file};
        if ($alignment_component ne "default") {
            $alignment_file_base .= "_$alignment_component";
        }
        my $alignment_file = $alignment_file_base . "_aligned.fa";
        $seq_al->write_fasta($alignment_file, $alignment_component);
        push @outputs, [$alignment_file, "aligned_${alphabet}_fasta"];
        if ($params->{recipe} eq 'phyml') {
            $alignment_file = $alignment_file_base . "_aligned.phy";
            $seq_al->write_phylip($alignment_file, $alignment_component);
            push @outputs, [$alignment_file, "phylip"];
        }
        print(STDERR "alignment file:  $alignment_file, size=" . -s $alignment_file . "\n") if $debug;

        my $tree_builder = new Tree_Builder($alignment_file);
        
        $tree_builder->set_program($params->{recipe});

        if ($model) {
           $tree_builder->set_model($model);
        }
        $tree_builder->set_output_base($alignment_file_base);
        #if ($params->{bootstrap}) {
        #    $tree_builder->set_bootstap_reps($params->{bootstrap});
        #}
        my $treeFile = $tree_builder->build_tree();
        $step_info->{command_line} .= "\n" if $step_info->{command_line};
        $step_info->{command_line} .= $tree_builder->get_analysis_commandline();
        $step_info->{details} .= "\n" if $step_info->{details};
        $step_info->{details} .= $tree_builder->get_analysis_stderrout();
        #my $logFile = $tree_builder->get_analysis_stderrout();
        push @outputs, ([$treeFile, "nwk"]);
    }
    end_step();

    #push @outputs, ([$logFile, "txt"]);
    # optionally re-write newick file with original sequence IDs if any were altered
    
    # generate tree graphic using figtree for all trees generated
    for my $file_record (@outputs) {
        if ($file_record->[1] eq 'nwk') {
            my $treeFile = $file_record->[0];
            print STDERR "About to call generate_tree_graphic($treeFile, $num_seqs, 'SVG')\n";
            my $tree_graphic = generate_tree_graphic($treeFile, $num_seqs, 'SVG');
            push @outputs, [$tree_graphic, 'SVG'];
            print STDERR "tree_file $treeFile\n";
        }
    }

    my %db_link_count;
    my $database_link = undef;
    my %data_source_count;
    my @ids = $seq_al->get_ids();
    for my $id (@ids) {
        my $link = $seq_al->get_metadata($id, 'database_link');
        if ($link) {
            $database_link = $link unless $database_link;
            $db_link_count{$link}++;
            if ($db_link_count{$link} > $db_link_count{$database_link}) {
                $database_link = $link;
            }
            my $data_source = $seq_al->get_metadata($id, 'data_source');
            if ($data_source) {
                $data_source_count{$data_source}++
            }
        }
    }
    
    my @command = ('p3x-newick-to-phyloxml', '-r', '[^(,)]+\_\@\_');
    if ($database_link) { # activate metadata retrieval from database
        push @command, ('-l', $database_link, '-g', join(',',@genome_metadata_fields), '-f', join(',', @feature_metadata_fields));
    }
    if ((scalar keys %data_source_count) > 1) {
        # write data sources to a tsv file and invoke adding it to phyloxml
        open F, ">data_source.tsv";
        print F "seq_id\tGroup\n";
        for my $id (@ids) {
            my $data_source = $seq_al->get_metadata($id, 'data_source'); 
            $data_source = "NA" unless $data_source;
            print F "$id\t$data_source\n";
        }
        close F;
        push @command, ("--annotationtsv", "data_source.tsv");
    }
    my ($step_comments, $step_info) = start_step("Format Tree to PhyloXML");
    for my $file_record (@outputs) {
        if ($file_record->[1] eq 'nwk') {
            my $treeFile = $file_record->[0];
            print STDERR "About to call p3x-newick-to-phyloxml on $treeFile\n" if $debug;
            print STDERR "run: " . join(' ', (@command, $treeFile)), "\n";
            $step_info->{command_line} = join(' ', (@command, $treeFile)) . "\n";
            my $ok = IPC::Run::run([@command, $treeFile]);
            my $phyloxml_file = $treeFile;
            $phyloxml_file =~ s/.nwk//;
            $phyloxml_file .= ".phyloxml";
            push @outputs, [$phyloxml_file, "phyloxml"];
        }
    }
    end_step();
     
    my $html_file = "$params->{output_file}_gene_tree_report.html";
    my $report_title = "Gene Tree Report";
    if ($params->{tree_type} eq 'viral_genome') {
        $html_file = "$params->{output_file}_virus_genome_tree_report.html";
        $report_title = "Virus Genome Tree Report";
    }
    my @tree_graphic_files;
    for my $file_record (@outputs) {
        if ($file_record->[1] eq 'SVG') {
            push @tree_graphic_files, $file_record->[0];
        }
    }
    write_report($html_file, $report_title, \@tree_graphic_files);
    push @outputs, [$html_file, "html"];

    print STDERR '\@outputs = '. Dumper(\@outputs);
    my $output_folder = $app->result_folder();
    for my $output (@outputs) {
        my($ofile, $type) = @$output;
        next if $type eq 'folder';
        
        if (! -f $ofile) {
            warn "Output file '$ofile' of type '$type' does not exist\n";
            next;
        }
        
        my $filename = basename($ofile);
        #print STDERR "Output folder = $output_folder\n";
        print STDERR "Saving $filename => $output_folder as $type\n" if $debug;
        if (0) { # for some reason this doesn't work
           $app->workspace->save_file_to_file($ofile, {}, "$output_folder/$filename", $type, 1,
               (-s $ofile > $shock_cutoff ? 1 : 0), # use shock for larger files
               $global_token);
        }
        else { # fall back to calling CLI
            my $ext = $1 if $ofile =~ /.*\.(\S+)$/;
            my @cmd = ("p3-cp", "-f", "-m", "${ext}=$type", $ofile, "ws:" . $app->result_folder);
            print STDERR "@cmd\n";
            my $ok = IPC::Run::run(\@cmd);
            if (!$ok)
            {
                warn "Error $? copying output with @cmd\n";
            }
        }
    }
    print STDERR "$tmpdir\n" if $debug;
    chdir($original_wd); # change back to the starting working directory
    my $time2 = `date`;
    print STDERR "Start: $time1\tEnd:   $time2\n";
    write_output("Start: $time1\tEnd:   $time2", "$tmpdir/DONE");
}

sub select_sequence_identifier {
    # fixed prioritization of elements to favor for best sequence identifier
    my ($seq_data) = @_;
    return $seq_data->{adjusted_id} if exists $seq_data->{adjusted_id};
    return $seq_data->{sequence_id} if exists $seq_data->{sequence_id};
    return $seq_data->{feature_id} if exists $seq_data->{feature_id};
    return $seq_data->{patric_id} if exists $seq_data->{patric_id};
    return $seq_data->{genome_id} if exists $seq_data->{genome_id};
    return $seq_data->{user_identifier} if exists $seq_data->{user_identifier};
    die "Couldn't find a usable sequence identifier for seq_data $seq_data";
}

sub retrieve_feature_metadata_by_patric_id {
    my ($api, $patric_id, $feature_fields) = @_;
    my $select_string = "select(" . join(",", @$feature_fields) . ")"; 
    $patric_id = uri_escape($patric_id);
    my $query = "eq(patric_id,$patric_id)&$select_string";
    my ($resp, $data) = $api->submit_query('genome_feature', $query);
    #print STDERR "response=$resp\n" if $debug;
    return $data->[0]; #only one element in return, which is a hash reference
}

sub retrieve_feature_metadata_by_feature_id {
    my ($api, $feature_id, $feature_fields) = @_;
    my $select_string = "select(" . join(",", @$feature_fields) . ")"; 
    #my $url = "$data_url/genome_feature/?eq(feature_id,($feature_id))&$select_string";
    #print STDERR "query=$url\n";
    #my $resp = curl_json($url);
    $feature_id = uri_escape($feature_id);
    my $query = "eq(feature_id,$feature_id)&$select_string";
    my ($resp, $data) = $api->submit_query('genome_feature', $query);
    #print STDERR "response=$resp\n" if $debug;
    return $data->[0]; #only one element in return, which is a hash reference
}

sub get_genome_metadata {
    my ($genome_ids, $fields) = @_;
    print STDERR "in get_genome_metadata: genome ids = ", join(", ", @$genome_ids), "\n" if $debug;
    my %genome_metadata = ();
    #return \() unless scalar(@$genome_ids);
    $fields = ['species'] unless $fields;
    my $select_string = "select(" . join(",", 'genome_id', @$fields) . ")"; 
    my @genome_ids_copy = @$genome_ids; # copy ids so we don't delete them with splice
    while (@genome_ids_copy) {
        my @id_sample = splice(@genome_ids_copy, 0, 20);
        my $escaped_ids = join(",", @id_sample);
        my $url = "$data_url/genome/?in(genome_id,($escaped_ids))&$select_string";
        print STDERR "query=$url\n";
        my $resp = curl_json($url);
        for my $member (@$resp) {
            #print STDERR "member: ", join(", ", sort(keys %$member)), "\n";
            my $genome_id = $member->{genome_id};
            for my $field (@$fields) {
                $genome_metadata{$genome_id}{$field} = $member->{$field} unless $field eq 'genome_id';
            }
        }
    }
    return \%genome_metadata
}

sub label_tree_with_metadata {
    my ($input_newick, $metadata, $label_fields) = @_;
    my ($step_comments, $step_info) = start_step("Label Tree With Metadata");
    my $comment = "label_tree_with_metadata called for fields: " . join(", ", @$label_fields);
    for my $field (@$label_fields) {
        print STDERR "metadata does not include $field\n" unless exists $metadata->{$field};
    }
    print STDERR $comment, "\n";
    my $output_newick = $input_newick;
    $output_newick =~ s/\..{2,6}$//; #remove extension
    $output_newick .= "_relabeled_" . join("_", @$label_fields) . ".nwk";
    my $newick = read_file($input_newick);
    print STDERR "Got newick string: ", substr($newick, 0, 50), "\n";
    my $start_field = $label_fields->[0];
    for my $seq_id (keys %{$metadata->{$start_field}}) {
        my @values = ();
        for my $field (@$label_fields) {
            push @values, $metadata->{$field}{$seq_id};
        }
        my $sub = join("|", @values);
        $seq_id =~ s/\|/\\\|/g;
        $sub =~ tr/[](),://d; #remove characters that break newick structure
        $sub =~ tr/ /_/; #remove characters that break newick structure
        print STDERR "substitute $seq_id with $sub\n";
        $newick =~ s/$seq_id/$sub/;
    }
    print STDERR "Modified newick string: ", substr($newick, 0, 50), "\n";
    write_file($output_newick, $newick);
    end_step("Label Tree With Metadata");
    return $output_newick;
}

sub generate_tree_graphic {
    my ($input_newick, $num_tips, $graphic_format) = @_;
    my ($step_comments, $step_info) = start_step("Generate Tree Graphic");
    my $file_base = basename($input_newick);
    $file_base =~ s/\..{2,6}//;
    my $tree_graphic_file = "$file_base." . lc($graphic_format);
    my $nexus_file = "$file_base.nex";
    my $comment = "run figtree input = $input_newick, output = $tree_graphic_file";
    #push @{$step_comments}, $comment;
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

    my @cmd = ("figtree", "-graphic", $graphic_format);
    
    if ($num_tips > 3) {
        my $height = 20 + 18 * ($num_tips); # this is an empirical correction factor to avoid taxon name overlap
        push @cmd, '-height', $height;
    }
    push @cmd, $nexus_file, $tree_graphic_file;
    $comment = join(" ", @cmd);
    $step_info->{command_line} = $comment;
    print STDERR "$comment\n";

    my ($stdout, $stderr) =  run_cmd(\@cmd);
    $step_info->{stdout} = $stdout;
    $step_info->{stderr} = $stderr;
    end_step();
    return $tree_graphic_file;
}

sub curl_text {
    my ($url) = @_;
    my @headers = ("-H", "Accept:text/tsv");
    if ($global_token) { push @headers, "-H", "Authorization:$global_token"}
    my @cmd = ("curl", $url, @headers);
    print STDERR join(" ", @cmd), "\n" if $debug;
    my ($out) = run_cmd(\@cmd);
    return $out;
}

sub curl_json {
    my ($url) = @_;
    my @headers = ("-H", "Accept:application/json");
    if ($global_token) { push @headers, "-H", "Authorization:$global_token"}
    my @cmd = ("curl", $url, @headers);
    print STDERR join(" ", @cmd), "\n" if $debug;
    my ($out) = run_cmd(\@cmd);
    print STDERR $out, "\n" if $debug;
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
    my ($id) = @_;
    # return $id; # testing
    my $ws = get_ws();
    my $token = get_token();
    
    my $base = basename($id);
    my $file = $base;
    # return $file; # testing
    print STDERR "get_ws_file:  base=$base, file=$file, ws=$ws, token=$token\n";
    
    my $fh;
    open($fh, ">", $file) or die "Cannot open $file for writing: $!";

    print STDERR "GET WS => $base $id\n";
    system("ls -lrta ");

    eval {
	$ws->copy_files_to_handles(1, $token, [[$id, $fh]]);
    };
    if ($@)
    {
        die "ERROR getting file $id\n$@\n";
    }
    close($fh);
    print "$id $file:\n";
    system("ls -lrta ");

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
