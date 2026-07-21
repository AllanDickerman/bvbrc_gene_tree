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
my @original_sequence_ids; # list of items requested, can be different from those actually obtained
my $sequence_identifier_type; # feature_id or genome_id or user_specified

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
    unless ($name eq $current_step) {
        print STDERR "Problem! at end_step name is wrong: $name should be $current_step";
    }
    my $step_index = scalar @analysis_step - 1;
    my $step_info = $analysis_step[$step_index];
    $step_info->{end_time} = time();
}

sub add_analysis_step { #allow adding a step recorded by Tree_Builder object
    my $step_info = shift;
    push @analysis_step, $step_info;
}

sub write_report {
    my ($output_file, $title, $tree_graphic_file) = @_;
    print STDERR "write_report()\n";
    open F, ">$output_file";
    print F "<HTML>\n<h1>$title</h1>\n";
    if (-e $tree_graphic_file) {
        my $element_name = 'tree_plot';
        print F "<script>\nfunction toggle_$element_name() {
          var x = document.getElementById(\"$element_name\");
            if (x.style.display == \"none\") {
                  x.style.display = \"block\";
              } else {
                  x.style.display = \"none\";
              }
          }\n</script>\n";
        print STDERR "<script>\nfunction toggle_$element_name() \n" if $debug;
        print F "FigTree Plot: <button onclick=\"toggle_$element_name()\">Show/Hide</button>\n";
        print F "<div id=\"$element_name\" style=\"display:block; background:#ffffff\" \n";
        print F "    onclick=\"toggle_$element_name()\">\n";
        my $svg_text = read_file($tree_graphic_file);
        print F $svg_text, "\n\n";
        print F "</div>\n";
    }
    else { print STDERR "Tree graphic file not found: $tree_graphic_file\n"; }
    print F "<h2>Analysis Steps</h2>\n";
    my $start_time = $analysis_step[0]->{start_time};
    my $time_string = localtime($start_time);
    my $duration = time() - $start_time;
    print F "<p>Start time $time_string<br>\n";
    print F "Duration: $duration\n";
    for my $step (@analysis_step) {
        print F "<h3>$step->{name}</h3>\n";
        if (exists $step->{command_line}) {
            print F "<B>command line:</b><br>\n";
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
        if (exists $step->{stdout} and length($step->{stdout}) > 5) {
            #accommodate steps from Tree_Builder
            $step->{details} = $step->{stdout};
        }
        if (exists $step->{details} and $step->{details} =~ /\S/) {
            my $element_name = "$step->{name}_details";
            $element_name =~ tr/ /_/;
            print F "<script>\nfunction toggle_$element_name() {
              var x = document.getElementById(\"$element_name\");
                if (x.style.display == \"none\") {
                      x.style.display = \"block\";
                  } else {
                      x.style.display = \"none\";
                  }
              }\n</script>\n";
            print STDERR "<script>\nfunction toggle_$element_name() \n" if $debug;
            print F "Details: <button onclick=\"toggle_$element_name()\">Show/Hide</button>\n";
            print F "<div id=\"$element_name\" style=\"display:none; background:#f0f0f0\" \n";
            print F "    onclick=\"toggle_$element_name()\">\n";
            print F "<pre>\n", $step->{details}, "\n</pre></div>\n";
        }
        my $duration = $step->{end_time} - $step->{start_time};
        if ($duration > 10) { 
            print F "<p>Duration ", $duration, " seconds.\n";
        }
    }
    print F "</HTML>\n";
}

sub retrieve_sequence_data {
    my ($app, $params, $api) = @_;
    my ($step_comments, $step_info) = start_step("Gather Sequence Data");
    my @master_seq_list; # concatenate seq_list for each data source, return this
    my %master_seq_ids; # control that sequence IDs are unique (modify user fast IDS, drop repeat feature or genomes)
    my $comment;
    my ($aligned_state, $any_in_memory) = (0, 0);
    $aligned_state = (scalar(@{$params->{sequences}}) == 1 and $params->{sequences}->[0]->{type} =~ /Aligned/i); 
    print STDERR "Number of sequence data sets = ", scalar(@{$params->{sequences}}), "\n";
    my $total_seqs = 0;
    my @empty_sequences; # keep track of entries lacking sequence data
    for my $sequence_source (@{$params->{sequences}}) {
        print STDERR "data item: $sequence_source->{type}, $sequence_source->{filename}\n";
        push @{$step_comments}, "reading $sequence_source->{type} $sequence_source->{filename}";
        my $num_seqs = 0;
        if ($sequence_source->{type} =~ /FASTA/i) {
            # then it is one of the fasta formats in the workspace
            my $local_file = $sequence_source->{filename};
            $local_file =~ s/.*\///; #remove path, leave filename
            print STDERR "About to copy file $sequence_source->{filename} to $local_file\n";
            $app->workspace->download_file($sequence_source->{filename}, $local_file, 1, $global_token);
            my @seq_list;
            open F, $local_file;
            my $seq_item;
            while (<F>) {
                if (/^>(\S+)/) {
                    my $user_identifier = $1;
                    push @original_sequence_ids, $user_identifier;
                    my %temp = ();
                    $seq_item = \%temp;
                    push @seq_list, $seq_item;
                    $seq_item->{id} = $user_identifier;
                    $seq_item->{user_identifier} = $user_identifier;
                    $seq_item->{data_source} = $local_file;
                    $seq_item->{sequence} = '';
                    if ($user_identifier =~ tr/[]():',/_______/) {
                        print STDERR "replacing identifier $seq_item->{id} with $user_identifier\n";
                        $seq_item->{id} = $user_identifier;
                    }
                    if (exists $master_seq_ids{$seq_item->{id}}) {
                        while (exists $master_seq_ids{$seq_item->{id}}) {
                            $seq_item->{id} .= "_d" # make unique by appending tag
                        }
                    }
                    $master_seq_ids{$seq_item->{id}} = 1; #remember this identifier and protect if from collisions
                    if ($user_identifier =~ /^fig\|\d+\.\d+\..{3}\.\d+$/) {
                        print "try user identifier as patric_id: $user_identifier\n" if $debug;
                        $seq_item->{database_link} = 'patric_id';
                    }
                    elsif ($user_identifier =~ /^\w+\.\d+\.\d+\.\w+\.\w+\.\d+\.\d+\.(fwd|rev)$/) {
                        print "try user identifier as feature_id: $user_identifier\n" if $debug;
                        $seq_item->{database_link} = 'feature_id';
                    }
                    elsif ($user_identifier =~ /(\d+\.\d+$)/) {
                        print "try user identifier as genome_id $user_identifier\n" if $debug;
                        $seq_item->{database_link} = 'genome_id';
                    }
                    elsif ($user_identifier =~ /(acc\|):?(\d+\.\d+$)/) {
                        print "try user identifier as genome_id $user_identifier\n" if $debug;
                        $seq_item->{database_link} = 'genome_id';
                    }
                    $num_seqs++;
                }
                else {
                    chomp;
                    $seq_item->{sequence} .= $_;
                }
            }
            close F;
            # save non-empty sequences to master seq list
            my $num_empty = 0;
            for my $item (@seq_list) {
                if (length($item->{sequence}) > 0) {
                    push @master_seq_list, $item;
                }
                else {
                    push @empty_sequences, $item->{id};
                    $num_empty++;
                }
            }
            $comment = "number of sequence entries: $num_seqs";
            push @{$step_comments}, $comment;
            print STDERR $comment, "\n";
            if ($num_empty) {
                $comment = "number lacking sequence data: $num_empty";
                push @{$step_comments}, $comment;
                print STDERR $comment, "\n";
            }
        }
        elsif ($sequence_source->{type} eq "feature_group") {
            # need to get feature sequences from database 
            my $feature_group = $sequence_source->{filename};
            $comment = "retrieving sequences for feature group $feature_group";
            print STDERR $comment, "\n";
            push @{$step_comments}, $comment;
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
            my $num_empty = 0;
            $feature_group =~ s/.*\///; # remove path preceding name of feature group
            for my $item (@$seq_list) {
                if (exists $master_seq_ids{$item->{feature_id}}) { # block repeats of same feature
                    print STDERR "got repeat of feature $item->{feature_id}, skipping.\n";
                }
                else {
                    push @original_sequence_ids, $item->{feature_id};
                    $item->{id} = $item->{feature_id};
                    $item->{data_source} = $feature_group;
                    $item->{database_link} = 'feature_id';
                    $master_seq_ids{$item->{feature_id}} = 1;
                    if (length($item->{sequence}) > 0) {
                        push @master_seq_list, $item;# save non-empty sequences to master seq list
                    }
                    else {
                        push @empty_sequences, $item->{id};
                        $num_empty++;
                    }
                }
            }
            $num_seqs = scalar @$seq_list;
            $comment = "number of sequence fetures retrieved: $num_seqs";
            push @{$step_comments}, $comment;
            print STDERR $comment, "\n";
            if ($num_empty) {
                $comment = "number lacking sequence data: $num_empty";
                push @{$step_comments}, $comment;
                print STDERR $comment, "\n";
            }
        }
        elsif ($sequence_source->{type} eq "genome_group") {
            if (0 and scalar @{$params->{sequences}} > 1) {
                print STDERR "Genome group $sequence_source->{filename} combined with other sequence inputs. This case is not handled yet. Exiting.\n";
                `cd ..`; # to allow temp directory to be deleted
                exit(1);
            }
            my $genome_group = $sequence_source->{filename};
            $comment = "retrieving sequences for genome group $genome_group\n";
            print STDERR "$comment\n";
            push @{$step_comments}, $comment;
            my $temp_genome_ids = $api->retrieve_patric_ids_from_genome_group($genome_group);
            my @genome_ids; # limit to non-redundant entries
            for my $id (@$temp_genome_ids) {
                if (exists $master_seq_ids{$id}) {
                    print STDERR " duplicate genome ID: $id, skipping.\n";
                }
                else {
                    print STDERR " genome $id\n";
                    push @genome_ids, $id;
                }

            }
            my @genome_validation_fields, ('genome_id', 'contigs', 'superkingdom', 'genome_length');
            my @genome_validation_data = $api->retrieve_genome_metadata(\@genome_ids, \@genome_validation_fields);
            print STDERR "examine genome metadata to test for single-sequence virus under $max_genome_length:\n";
            print STDERR join("\t", @genome_validation_fields), "\n";
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
                    print STDERR "Problem: genome $info->{'genome_id'} has multiple contigs ($info->{'contigs'}).\nExiting.\n";
                    `cd ..`;
                    exit(1);
                }
            }
            print STDERR "All genomes are viruses, all have a single sequence, all are under $max_genome_length bases.\n";
            push @original_sequence_ids, @genome_ids;
            $num_seqs = scalar @genome_ids;
            my $num_empty = 0;
            $genome_group =~ s/.*\///; # remove path preceding name of genome group
            for my $item (@genome_validation_data) {
                my $genome_id = $item->{genome_id};
                my ($resp, $data) = $api->submit_query('genome_sequence', "eq(genome_id,$genome_id)", "sequence");
                #print "for $genome_id: resp = $resp\tdata=$data\tdata->[0]=$data->[0]\n" if $debug;
                $item->{sequence} = $data->[0]->{sequence};
                $item->{id} = $genome_id;
                $item->{data_source} = $genome_group;
                $item->{database_link} = 'genome_id';
            # save non-empty sequences to master seq list
                if (length($item->{sequence}) > 0) {
                    push @master_seq_list, $item;
                }
                else {
                    push @empty_sequences, $item->{id};
                    $num_empty++;
                }
            }
            $comment = "number of genomes retrieved: $num_seqs";
            push @{$step_comments}, $comment;
            print STDERR $comment, "\n";
            if ($num_empty) {
                $comment = "number lacking sequence data: $num_empty";
                push @{$step_comments}, $comment;
                print STDERR $comment, "\n";
            }

        }
        elsif ($sequence_source->{type} eq "feature_ids") {
            # need to get feature sequences from database 
            $comment = "retrieving sequences for feature ids";
            $comment .= ", this functionality in testing";
            push @{$step_comments}, $comment;
            my @feature_ids;
            for my $feature_id (@{$sequence_source->{feature_ids}}) {
                if (exists $master_seq_ids{$feature_id}) { # block repeats of same feature
                    print STDERR "got repeat of feature $feature_id, skipping.\n";
                }
                else {
                    push @feature_ids, $feature_id;
                    $master_seq_ids{$feature_id} = 1;
                }
            }
            print STDERR "\tfeature_ids = ", join(", ", @feature_ids), "\n" if $debug;
            my $query="in(feature_id,(" . join(',', @feature_ids) . "))";
            my ($req, $seq_list) = $api->submit_query('genome_feature', $query);
            my @md5_list;
            my $md5_type = ('aa_sequence_md5', 'na_sequence_md5')[$params->{alphabet} eq 'DNA'];
            for my $item  (@$seq_list) {
                push @md5_list, $item->{$md5_type};
            }
            my $seqs = $api->lookup_sequence_data_hash(\@md5_list);
            my $num_empty = 0;
            for my $item  (@$seq_list) {
                $item->{sequence} = $seqs->{$item->{$md5_type}};
                $item->{"database_link"}="feature_id"; 
                if (length($item->{sequence}) > 0) {
                    push @master_seq_list, $item;
                }
                else {
                    push @empty_sequences, $item->{id};
                    $num_empty++;
                }
            }
            $num_seqs = scalar keys %{$sequence_source->{sequences}};
            $comment = "number of features retrieved: $num_seqs";
            push @{$step_comments}, $comment;
            print STDERR $comment, "\n";
            if ($num_empty) {
                $comment = "number lacking sequence data: $num_empty";
                push @{$step_comments}, $comment;
                print STDERR $comment, "\n";
            }
        }
    }
    if (scalar @empty_sequences) {
        $comment = "records lacking sequence data: " . join(", ", @empty_sequences);
        push @{$step_comments}, $comment;
        print STDERR "$comment\n";
    }
    
    my $num_seqs = scalar @master_seq_list;
    print STDERR "near end of retrieve_sequence_data, number of sequences is $num_seqs\n" if $debug;

    $comment = $aligned_state ? "sequences are aligned" : "sequences need aligning";
    push @{$step_comments}, $comment;
    print STDERR "$comment\n";
    end_step("Gather Sequence Data");
    return \@master_seq_list;
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
    chdir($tmpdir); # do all work in temporary directory
   
    my @feature_metadata_fields = @default_feature_metadata_fields;
    if (exists $params->{feature_metadata_fields}) {
        @feature_metadata_fields = @{$params->{feature_metadata_fields}};
    }
    #ensure that feature_id and genome_id are retrieved
    push @feature_metadata_fields, "feature_id" unless grep(/feature_id/, @feature_metadata_fields);
    push @feature_metadata_fields, "genome_id" unless grep(/genome_id/, @feature_metadata_fields);

    my @genome_metadata_fields = @default_genome_metadata_fields;
    if (exists $params->{genome_metadata_fields}) {
        @genome_metadata_fields = @{$params->{genome_metadata_fields}};
    }

    my $seq_list = retrieve_sequence_data($app, $params, $api);
    my $seqids_are_genome_ids = 0; # indicate whether seq identifiers are links to BVBRC genomes in database
    my $database_link_type = undef;
    my $num_seqs = scalar @$seq_list;
    print STDERR "After retrieval, number of sequences is $num_seqs\n" if $debug;
    if ($num_seqs < 4) { #need at least 4 seuqences to build a tree
        print STDERR "After retrieval, number of sequences is $num_seqs, less than 4. Cannot build a tree.\n";
        exit(1);
    }

    my $is_aligned = 0;
    my $unaligned_fasta_file = "$params->{output_file}_unaligned.fa";
    my $aligned_fasta_file = "$params->{output_file}_aligned.fa";
    my $outfile;
    if (scalar @{$params->{sequences}} == 1 and $params->{sequences}->[0]->{type} =~ /Aligned/i) {
        open $outfile, ">$aligned_fasta_file";
        $is_aligned = 1;
    }
    else {
        open $outfile, ">$unaligned_fasta_file" or die "could not open fasta file for output";
    }
    my $database_link = undef;
    my %seqs_written;
    for my $seq_item (@$seq_list) {
        my $seq_id = $seq_item->{id};
        #debugging
        warn "multiple occurrence of $seq_id\n%{$seq_item}\n" if $seqs_written{$seq_id};
        $seqs_written{$seq_id} = 1;

        print $outfile ">$seq_id\n";
        #print STDERR "writing sequence for $seq_id\n";
        my $sequence = $seq_item->{sequence};
        $sequence =~ tr/-//d unless $is_aligned;
        print $outfile $sequence, "\n";
        if (exists $seq_item->{database_link}) {
            if ($database_link and $database_link ne $seq_item->{database_link}) {
                print STDERR "Two different database_links available: $database_link ne $seq_item->{database_link}\n";
            }
            $database_link = $seq_item->{database_link};
        }
    }
    close $outfile;
    my $alignment_modified = 0; # flag whether alignment is different from input data 
    unless ($is_aligned) {
        #run_muscle($unaligned_file, $aligned_fasta_file);
        run_mafft($unaligned_fasta_file, $aligned_fasta_file);
        $alignment_modified = 1;
    }
    if ($params->{trim_threshold} or $params->{gap_threshold})
    {
        my $trimmed_aligned_file = $aligned_fasta_file;
        $trimmed_aligned_file =~ s/.afa//;
        $trimmed_aligned_file =~ s/.fa//;
        $trimmed_aligned_file =~ s/.fasta//;
        $trimmed_aligned_file .= "_trimmed.afa";
        my $retval = trim_alignment($aligned_fasta_file, $trimmed_aligned_file, $params->{trim_threshold}, $params->{gap_threshold});
        if ($retval) {
            $alignment_modified = 1;
            print STDERR "trimmed aligned file written to $trimmed_aligned_file\n" if $debug;
            $aligned_fasta_file = $trimmed_aligned_file;
        }
    }
    my $alphabet = $params->{alphabet};
    push @outputs, [$aligned_fasta_file, "aligned_${alphabet}_fasta"] if $alignment_modified;
    run("echo $tmpdir && ls -ltr $tmpdir") if $debug;

    my $model = "LG"; # default for protein
    if ($params->{substitution_model}) {
        $model = $params->{substitution_model}
    }
    elsif ($params->{alphabet} =~ /DNA/i) {
        $model = "GTR"
    }
    my $recipe = "raxml"; #default
    if (defined $params->{recipe} and $params->{recipe}) {
        $recipe = lc($params->{recipe})
    }
    print STDERR "About to call tree program $recipe\n";
    my @tree_outputs;

    my $threads = 2;
    if (exists $ENV{P3_ALLOCATED_CPU}) {
        $threads = $ENV{P3_ALLOCATED_CPU};
        print STDERR "P3_ALLOCATED_CPU = $ENV{P3_ALLOCATED_CPU}\n";
    }
    my $tree_builder = new Tree_Builder($aligned_fasta_file, $alphabet, $threads);

    if ($model) {
       $tree_builder->set_model($model);
    }
    $tree_builder->set_output_base($params->{output_file}) if defined $params->{output_file};
    my $treeFile;
    my $bootstrap = $params->{bootstrap};
    if ($recipe eq 'raxml') {
        $treeFile = $tree_builder->build_raxml_tree($bootstrap);
    } elsif ($recipe eq 'phyml') {
        $treeFile = $tree_builder->build_phyml_tree($bootstrap);
    } elsif ($recipe eq 'fasttree') {
        $treeFile = $tree_builder->build_fasttree($bootstrap);
    } else {
        die "Unrecognized program: $recipe \n";
    }
    for my $index (0..$tree_builder->get_num_analysis_steps()-1) {
        my $step = $tree_builder->get_analysis_step($index);
        push @analysis_step, $step;
    }
    my $logFile = $tree_builder->get_log_file();
    push @outputs, ([$treeFile, "nwk"], [$logFile, "txt"]);
    # optionally re-write newick file with original sequence IDs if any were altered
    # not sure we want to do that as it could screw up downstream use of newick data i
    # (ideally, bad characters should be escaped, but does raxml know about escaped characters?)
    # generate tree graphic using figtree
    my $graphic_format = 'SVG';
    print STDERR "About to call generate_tree_graphic($treeFile, $num_seqs, $graphic_format)\n";
    my $tree_graphic = generate_tree_graphic($treeFile, $num_seqs, $graphic_format);
    push @outputs, [$tree_graphic, $graphic_format];
    print STDERR "tree_file $treeFile\n";
    if ($database_link) { # use system call to p3x-newick-to-phyloxml 
        my @command = ('p3x-newick-to-phyloxml', '-r', '[^(,)]+\_\@\_', '-l', $database_link, '-g', join(',',@genome_metadata_fields), '-f', join(',', @feature_metadata_fields)); # minus tree file (add later)
        my %data_source_count;
        for my $seq_item (@$seq_list) {
            my $data_source = $seq_item->{data_source};
            if ($data_source) {
                $data_source_count{$data_source}++
            }
        }
        if ((scalar keys %data_source_count) > 1) {
            # write data sources to a tsv file and invoke adding it to phyloxml
            open F, ">data_source.tsv";
            print F "seq_id\tGroup\n";
            for my $seq_item (@$seq_list) {
                my $data_source = "NA";
                if (exists $seq_item->{data_source}) {
                    $data_source = $seq_item->{data_source};
                }
                print F "$seq_item->{id}\t$data_source\n";
            }
            close F;
            push @command, ("--annotationtsv", "data_source.tsv");
        }
        push @command, $treeFile;
            
        #debugging
        print STDERR "execute system call: (as array):\n" . join(' ', @command), "\n";
        system(@command);

        my $phyloxml_file = undef;
        opendir my $dir_handle, '.' or die "Couldn't open dir '.': $!";
        my @files = readdir $dir_handle;
        for my $file (@files) {
            $phyloxml_file = $file if ($file =~ /.phyloxml$/);
        }
        push @outputs, [$phyloxml_file, "phyloxml"];
    }
    
    my $html_file = "$params->{output_file}_gene_tree_report.html";
    my $report_title = "Gene Tree Report";
    if ($params->{tree_type} eq 'viral_genome') {
        $html_file = "$params->{output_file}_virus_genome_tree_report.html";
        $report_title = "Virus Genome Tree Report";
    }
    write_report($html_file, $report_title, $tree_graphic);
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

sub trim_alignment {
    my ($aligned_fasta_file, $trimmed_aligned_file, $trim_threshold, $gap_threshold) = @_;
    my ($trim_comments, $trim_info) = start_step("Trim Alignment");
    my $comment = "performing trimming on alignment: trim_threshod=$trim_threshold, gap_threhold=$gap_threshold";
    print STDERR "$comment\n";
    push @{$trim_comments}, $comment;
    my $alignment = new Sequence_Alignment($aligned_fasta_file);
    print STDERR "trim_alignment: alignment=$alignment\n" if $debug;
    $trim_info->{details} = "Before trimming:\n" .  $alignment->write_stats();
    my $alignment_changed = 0;
    if ($trim_threshold > 0) {
        $comment = "trim ends of alignment to density $trim_threshold";
        print STDERR "$comment\n";
        push @{$trim_comments}, $comment;
        my ($left_trim_cols, $right_trim_cols) = $alignment->end_trim($trim_threshold);
        if ($left_trim_cols or $right_trim_cols) {
            $comment = "ends trimmed: $left_trim_cols columns on left, $right_trim_cols columns on right.";
            print STDERR "$comment\n";
            push @{$trim_comments}, $comment;
            $alignment_changed = 1;
        }
        else {
            $comment = "no end-columns needed to be trimmed";
            push @{$trim_comments}, $comment;
            print STDERR "$comment\n";
        }
    }
    if ($gap_threshold > 0) {
        $comment = "delete any sequences with gap proportion greater than $gap_threshold";
        print STDERR "$comment\n";
        push @{$trim_comments}, $comment;
        my $deleted_seqs = $alignment->delete_gappy_seqs($gap_threshold);
        if (scalar @$deleted_seqs) {
            $comment = scalar(@$deleted_seqs) . " gappy sequences deleted: ". join(", ", @$deleted_seqs);
            print STDERR "$comment\n";
            push @{$trim_comments}, $comment;
            $alignment_changed = 1;
        }
        else {
            $comment = "no sequences needed to be deleted";
            push @{$trim_comments}, $comment;
            print STDERR "$comment\n";
        }
    }
    if ($alignment_changed) {
        $trim_info->{details} .= "After trimming:\n" . $alignment->write_stats();
    }
    my $total_seqs = $alignment->get_ntaxa();
    if ($total_seqs < 4) {
        $comment = "total sequences: $total_seqs, too few to build a tree (minimum 4)";
        push @{$trim_comments}, $comment;
        print STDERR "$comment\n";
    }
    if ($alignment_changed) {
        $alignment->write_fasta($trimmed_aligned_file); 
        $comment = "writing trimmed alignment to $trimmed_aligned_file\n";
        push @{$trim_comments}, $comment;
        print STDERR $comment;
    }
    end_step("Trim Alignment");
    return $alignment_changed;
}

sub run_raxml {
    my ($alignment_file, $alphabet, $model, $output_name, $bootstrap) = @_;
    my ($step_comments, $step_info) = start_step("Phylogenetic Inference with RAxML");
    print STDERR "In run_raxml (with RELL support), alignment = $alignment_file\n";
    my $parallel = $ENV{P3_ALLOCATED_CPU};
    $parallel = 2 if $parallel < 2;
    
    my $cwd = getcwd();
    $model = uc($model); 
    if ($alphabet eq 'DNA') {
        $model = 'GTRGAMMA'
    }
    else {
        $model = 'LG' if $model !~ /DAYHOFF|DCMUT|JTT|MTREV|WAG|RTREV|CPREV|VT|BLOSUM62|MTMAM|LG|MTART|MTZOA|PMB|HIVB|HIVW|JTTDCMUT|FLU|STMTREV|DUMMY|DUMMY2|AUTO|LG4M|LG4X|PROT_FILE|GTR_UNLINKED|GTR/i;
        $model = "PROTCAT". $model;
    }

    my @cmd = ("raxmlHPC-PTHREADS-SSE3");
    push @cmd, ("-T", $parallel);
    push @cmd, ("-p", "12345");
    push @cmd, ("-m", $model);
    push @cmd, ("-s", basename($alignment_file));
    push @cmd, ("-n", $output_name);
    my $support_trees;
    if ($bootstrap) {
        push @cmd, ("-f",  "d"); # just do ML search (generate bootstrap trees separately later)
        $support_trees = "RAxML_bootstrap." . $output_name . "_bootstrap";
    }
    else {
        push @cmd, ("-f", "D"); # generate RELL replicate trees
        $support_trees = "RAxML_rellBootstrap." . $output_name;
    }
    
    my $comment = "command = ". join(" ", @cmd);
    push @{$step_comments}, $comment;
    print STDERR "$comment\n\n";
   
    my ($out, $err) = run_cmd(\@cmd);
    print STDERR "STDOUT:\n$out\n";
    print STDERR "STDERR:\n$err\n";
    #$step_info->{details} = $out;

    if ($bootstrap) {
        if ($alphabet eq 'DNA') {
            $model = 'GTRCAT'
        }
        @cmd = ("raxmlHPC-PTHREADS-SSE3");
        push @cmd, ("-T", $parallel);
        push @cmd, ("-m", $model);
        push @cmd, ("-s", basename($alignment_file));
        push @cmd, ("-n", $output_name . "_bootstrap");
        push @cmd, ("-b",  "12345", "-#", $bootstrap, "-p", "12345");
        ($out, $err) = run_cmd(\@cmd);
    }
    # now map replicate support numbers onto ML tree (from either bootstrap or RELL output)
    @cmd = ("raxmlHPC-PTHREADS-SSE3", "-f", "b", "-m", "GTRCAT"); #to map RELL bootstrap support onto ML tree
    push @cmd, ("-t", "RAxML_bestTree.".$output_name);
    push @cmd, ("-z", $support_trees);
    push @cmd, ("-n", $output_name . "_support");
    print STDERR "Map support values onto ML tree:\n", join(" ", @cmd), "\n";
    my ($out, $err) = run_cmd(\@cmd);
    print STDERR "STDOUT:\n$out\n";
    print STDERR "STDERR:\n$err\n";
    #$step_info->{details} .= $out;
    
    my @outputs;
    my $treeFile = $output_name . "_raxml_tree.nwk";
    move("RAxML_bipartitions.".$output_name. "_support", $treeFile);
    my $logFile = $output_name . "_raxml_log.txt";
    move("RAxML_info.".$output_name . "_support", $logFile);
    #my $details = read_file($logFile);
    #$step_info->{details} .= $details;
    push @outputs, [$treeFile, 'nwk'];
    push @outputs, [$logFile, 'txt'];

    run("ls -ltr");
    #end_step("Phylogenetic Inference with RAxML");
    return @outputs;
}

sub run_phyml {
    my ($alignment_file, $alphabet, $model, $output_name, $bootstrap) = @_;
    #my ($step_comments, $step_info) = start_step("Phylogenetic Inference with Phyml");

    my $cwd = getcwd();
    
    my $datatype = 'aa';
    if ($alphabet =~ /DNA/i) {
        $datatype = 'nt';
        $model = 'GTR' if $model !~ /HKY85|JC69|K80|F81|F84|TN93|GTR/;
    }
    else {
        $model = 'LG' if $model !~ /WAG|JTT|MtREV|Dayhoff|DCMut|RtREV|CpREV|VT|AB|Blosum62|MtMam|MtArt|HIVw|HIVb/;
    }

    my @cmd = ("phyml");
    push @cmd, ("-i", $alignment_file);
    push @cmd, ("-d", $datatype);
    push @cmd, ("-m", $model);
    if ($bootstrap) {
        push @cmd, ("-b", $bootstrap); # normal bootstrap replicates
    }
    else {
        push @cmd, ("-b", '-3'); # -1 gives approximate Likelihood Ratio Tests
        # -2 give Chi2-based parametric branch supports
        # -3 gives Shimodaira-Hasegawa (SH) support values
    }
    
    my $comment = "command = ". join(" ", @cmd);
    #push @{$step_comments}, $comment;
    print STDERR "$comment\n\n";
   
    my ($out, $err) = run_cmd(\@cmd);
    print STDERR "STDOUT:\n$out\n";
    print STDERR "STDERR:\n$err\n";
    #$step_info->{details} = $out;
    #$step_info->{stderr} = $err;
    #my $comment = "return code = $rc\n";
    #push @{$step_comments}, $comment;
    
    my @outputs;
    my $treeFile = $output_name."_phyml_tree.nwk";
    move($alignment_file."_phyml_tree.txt", $treeFile);
    my $logFile = $output_name."_phyml_log.txt";
    move($alignment_file."_phyml_stats.txt", $logFile);
    my $details = read_file($logFile);
    #$step_info->{details} .= $details;
    push @outputs, [$treeFile, 'nwk'];
    push @outputs, [$logFile, 'txt'];

    run("ls -ltr");

    #end_step("Phylogenetic Inference with Phyml");
    return @outputs;
}

sub run_fasttree {
    my ($alignment_file, $alphabet, $model, $output_name, $bootstrap) = @_;
    #my ($step_comments, $step_info) = start_step("Phylogenetic Inference with FastTree");
    my $treeFile = $output_name."_fasttree.nwk";
    my @cmd = ("FastTree", "-out", $treeFile);
    if ($alphabet =~ /DNA/i) {
        push @cmd, "-nt", "-gtr";
    }
    else {
        $model = 'lg' if $model !~ /lg|wag/i;
        $model = lc($model);
        push @cmd, "-$model";
    }

    push @cmd, $alignment_file;
    my $comment = join(" ", @cmd);
    #push @{$step_comments}, $comment; 
    print STDERR $comment . "\n\n";
   
    my ($out, $err) = run_cmd(\@cmd);
    print STDERR "STDOUT:\n$out\n";
    print STDERR "STDERR:\n$err\n";

    if ($bootstrap) {
        # use raxml to generate 100 data matrices
        # use fasttree to analyze them
        # use compareToBootstrap to count the per-branch support (as proportion)
        @cmd = "raxmlHPC-PTHREADS-SSE3  -f j -b 123 -# $bootstrap -s $alignment_file -n boot_matrix -m GTRCAT";
        print STDERR "Generate bootstrap matrices:\n" . join(" ", @cmd) . "\n";
        system(@cmd);
        die "raxml failed to produce BS1" unless -f "$alignment_file.BS1";
        my $multi_alignment_file = $alignment_file;
        $multi_alignment_file =~ s/\.{3,7}$//;
        $multi_alignment_file .= "_${bootstrap}BS.phy";
        @cmd = ("cat", "$alignment_file.BS*", ">", $multi_alignment_file);
        print STDERR "Concatenate into one file:\n" . join(" ", @cmd) . "\n";
        system(join(" ", @cmd)); # requires shell interpolation
        system("rm $alignment_file.BS*"); #allow shell expansion
        my $multi_tree_file = "multiple_bootstrap_trees.nwk";
        @cmd = ("FastTree", "-n", $bootstrap, "-out", $multi_tree_file);
        if ($alphabet =~ /DNA/i) {
            push @cmd, "-nt", "-gtr";
        }
        else {
            push @cmd, "-$model";
        }
        push @cmd, $multi_alignment_file;
        print STDERR "Generate trees for each bootstrapped data matrix:\n" . join(" ", @cmd) . "\n";
        system(@cmd);
        die "FastTree did not generate multiple tree file from boostrapped matrices" unless -f $multi_tree_file;
        my $tree_with_support = $output_name."_fasttree_boostrap_prop.nwk"; 
        @cmd = ("CompareToBootstrap", $treeFile, $multi_tree_file, ">", $tree_with_support);
        print STDERR "Map support onto best tree:\n" . join(" ", @cmd) . "\n";
        system(join(' ', @cmd));
        die "CompareToBootstrap did not generate tree with support values" unless -s $tree_with_support;
        $treeFile = $tree_with_support;
    }

    #$step_info->{details} = $err;
    
    my @outputs;
    my $logFile = $output_name."_fasttree_log.txt";
    write_file($logFile, $err);
    push @outputs, [$treeFile, 'nwk'];
    push @outputs, [$logFile, 'txt'];

    run("ls -ltr");

    #end_step("Phylogenetic Inference with FastTree");
    return @outputs;
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
    
    if ($num_tips > 40) {
        my $height = 600 + 15 * ($num_tips - 40); # this is an empirical correction factor to avoid taxon name overlap
        push @cmd, '-height', $height;
    }
    push @cmd, $nexus_file, $tree_graphic_file;
    $comment = join(" ", @cmd);
    add_step_command_line(join(" ", @cmd));
    print STDERR "$comment\n";

    my ($stdout, $stderr) =  run_cmd(\@cmd);
    $step_info->{stdout} = $stdout;
    $step_info->{stderr} = $stderr;
    end_step("Generate Tree Graphic");
    return $tree_graphic_file;
}

sub run_muscle {
    my ($unaligned, $aligned) = @_;
    print STDERR "run_muscle($unaligned, $aligned)\n";
    my ($step_comments, $step_info) = start_step("Align with muscle");
    my $cmd = ["muscle", "-in", $unaligned, "-out", $aligned];
    my $comment = join(" ", @$cmd);
    add_step_command_line($comment);
    print STDERR $comment, "\n" if $debug;
    my ($stdout, $stderr) =  run_cmd($cmd);
    $step_info->{details} = $stderr;
    end_step("Align with muscle");
}

sub run_mafft {
    my ($unaligned, $aligned) = @_;
    print STDERR "run_mafft($unaligned, $aligned)\n";
    my ($step_comments, $step_info) = start_step("Align with mafft");
    my $parallel = $ENV{P3_ALLOCATED_CPU};
    my $cmd = ["mafft", "--auto"];
    if ($parallel) {
        push @$cmd, "--thread", $parallel;
    }
    push @$cmd, $unaligned;
    my $comment = join(" ", @$cmd);
    add_step_command_line($comment);
    print STDERR $comment, "\n" if $debug;
    my $out;
    open $out, ">$aligned";
    run($cmd, ">", $out);
    end_step("Align with mafft");
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
