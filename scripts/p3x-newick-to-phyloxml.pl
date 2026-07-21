#!/usr/bin/env perl
=head1 Generate PhyloXML version of a tree with metadata for tip nodes

    p3x-generate-phyloxml.pl [options] newickFile

This script reads a newick tree file and writes the tree in PhyloXML format.
Optionally, it can add annotation describing the tips of the tree.
Annotation can be conveyed in a TSV file or fetched via the BVBRC data api.

=head2 Parameters

The one positional parameter is the path of the newick file.

The command-line options are as follows.

=over 4

=item annotationtsv

A tab-delimited file where the first column contains the tip IDs on the tree.
The column headings are the names of the fields.
The values are annotations for each field (column) for each tip (row).
An optional second header row with the first field 'provenance' will be used to prefix the phy:Property fields in the PhyloXML structure (othersize the value of annotationProvenance is used).

=item databaselink 

Identifies how the tree tip identifiers are to be interpreted when linking to the BVBRC database.
Options are 'patric_id', 'genome_id', 'feature_id'.

=item featurefields

Comma separated list of feature-level fields requested for annotating tips.

=item genomefields

Comma separated list of genome-level fields requested for annotating tips.

=back

=cut

use strict;
use P3DataAPI;
use P3Utils;
use URI::Escape;
use Cwd;
use File::Temp;
use Phylo_Tree; # should be in lib directory

#$| = 1;
# Get the command-line options.

my($opt, $usage) = P3Utils::script_opts('newickFile',
                ['databaselink|link|l=s', 'Name of database field that tree identifiers map to (feature_id, genome_id, patric_id).'],
                ['genomefields|g=s', 'Comma-separated list of genome fields to annotate each tree tip.', {default => 'species,strain,geographic_group,isolation_country,host_group,host_common_name,collection_year,genus,mlst'}],
                ['featurefields|f=s', 'Comma-separated list of feature fields to annotate each tree tip.', {default => 'product,accession'}],
                ['output_name=s', 'Output filename.'],
                ['name_with_all_fields', 'Construct output file to include all fields.'],
                ['overwrite|o', 'Overwrite existing files if any.'],
                ['verbose|debug|v', 'Write status messages to STDERR.'],
                ['name=s', 'Name for tree in phyloxml.'],
                ['remove_substring|r=s', 'Substring to be removed from finalized data'],
                ['description=s', 'Description of tree in phyloxml.'],
                ['id_translation=s', 'Translation table for tip IDs.'],
                ['annotationtsv|a=s', 'TSV file containing annotation for tips on the tree (or comma-delimited list)'],
                ['provenance|p=s', 'Provenance of annotation. (default: BVBRC)', { default => 'BVBRC' }],
        );

# Check the parameters.

if ($opt->verbose) {
    print "args=", join(", ", @ARGV), "\n";
    print "opt = %$opt\n";
    for my $key (keys %$opt) {
        print "\t$key\t$opt->{$key}\n";
    }
}

die($usage->text) if @ARGV != 1;

my $newickFile = shift;

# Get the debug flag.
my $debug = $opt->verbose;
if ($debug) {
    $debug=3;
    Phylo_Tree::set_debug($debug); 
    Phylo_Node::set_debug($debug); 
    print STDERR "args = ", join("\n", @ARGV), "\n";
}

my $tmpdir = File::Temp->newdir( "/tmp/TreeAnnotation_XXXXX", CLEANUP => !$debug );
system("chmod", "755", "$tmpdir");
print STDERR "created temp dir: $tmpdir, cleanup = ", !$debug, "\n" if $debug;
my $original_wd = getcwd();
my $workspace_dir;
my @fields = split("/", $newickFile);
if (scalar @fields > 2 and $fields[1] =~ '@') {
    print STDERR "looking for $newickFile in user workspace\n" if $opt->verbose;
    # probably a user workspace path
    chdir($tmpdir); # do all work in temporary directory
    my ($user, $brc) = split('@', $fields[1]);
    if ($brc and $brc =~ /patricbrc\.org|bvbrc/) {
        my $workspace_newick = $newickFile;
        $newickFile = pop @fields; # grab the last '/'-delimited field in user workspace path
        $workspace_dir = join('/', @fields);
        print STDERR "Parsed user=$user, brc=$brc, file=$newickFile\n" if $debug;
        my $ls_result = `p3-ls '$workspace_newick'`;
        chomp $ls_result;
        if ($ls_result ne $workspace_newick) {
            print "'$workspace_newick'\n'$ls_result'\n" if $opt->verbose;
            print "Cannot access $workspace_newick\nPerhaps not logged in as user $user\n";
            exit(1);
        }
        if (-f $newickFile and not $opt->overwrite) {
            print "Refusing to overwrite local file $newickFile. Use --overwrite (or -f) to enable overwrite.\n";
            exit(1);
        }
        my @command = ('p3-cp');
        push @command, '-f' if $opt->overwrite;
        push @command, "ws:" . $workspace_newick;
        push @command, '.';
        print STDERR "command to copy tree from workspace:\n@command\n" if $opt->verbose;
        my $rc = system(@command);
	die "Failure $rc running @command" unless $rc == 0;
        unless (-f $newickFile) {
            die "Failed to copy $workspace_newick to local file system.";
        }
    }
}
unless (-f $newickFile) {
    print "Cannot find file $newickFile.";
    exit(1);
}

my $tree = new Phylo_Tree(); 
#print STDERR "read tree. Newick is\n", $tree->write_newick(), "\n" if $debug;

$tree->read_newick($newickFile);

if ($opt->databaselink) {
    $tree->set_type($opt->databaselink);
}
if ($opt->name) {
    $tree->set_name($opt->name);
}
if ($opt->description) {
    $tree->set_description($opt->description);
}
if ($opt->id_translation) { #translate tip IDs
    print STDERR "reading tip label translation from $opt->id_translation\n" if $debug;
    open F, $opt->id_translation;
    my %gtdb_bvbrc;
    while (<F>) {
        chomp;
        my ($gtdb, $bvbrc) = split("\t");
        $gtdb_bvbrc{$gtdb} = $bvbrc;
    }
    # iterate over tips of tree
    # swap GTDB ID for BVBRC genome ID
    # store GTDB ID as metadata for node
    print STDERR "swapping tip IDs\n" if $debug;
    $tree->swap_tip_names(\%gtdb_bvbrc, 'original_id');
}

my %meta_column; #first key is column (field name), second key is row (tip ID)
if ($opt->annotationtsv) {
    for my $metadata_file ( split(',', $opt->annotationtsv) ) {
        print STDERR "reading metadata from $metadata_file\n" if $debug;
        open F, $metadata_file;
        $_ = <F>;
        chomp;
        my @header = split("\t");
        shift @header; #remove first element
        for my $column_head (@header) {
            if (exists $meta_column{$column_head}) {
                print STDERR "metadata column $column_head from file $metadata_file already exists: data may be over-written\n";
            }
            $meta_column{$column_head}{'column_head'} = $column_head;
        }
        print STDERR "Header fields: " . join(", ", @header) . "\n" if $debug;
        while (<F>) {
            s/^#//; # remove leading pound sign, if any
            chomp;
            my @fields = split("\t");
            my $id = shift @fields; # remove first element
            print STDERR " metadata for $id:" if $debug;
            for my $i (0 .. $#header) {
                my $column_head = $header[$i];
                my $val = $fields[$i];
                $meta_column{$column_head}{$id} = $val;
                print STDERR "\t$column_head=$val" if $debug;
            }
            print STDERR "\n" if $debug;
        }
    }
    for my $field (keys %meta_column) {
        $tree->add_tip_annotation($field, $meta_column{$field});
    }
}

if ($opt->databaselink) {
    # Get access to PATRIC.
    my $api = P3DataAPI->new();
    my $treeIds = $tree->get_tip_names();
    my $num_tips = scalar @$treeIds;
    my $limit = "limit($num_tips)";
    my $link = $opt->databaselink;
    print STDERR "tree IDs are: ", join(", ", @$treeIds), "\n" if $debug;
    my @escaped_IDs = map { uri_escape $_ } @$treeIds;
    my $query = "in($link,(" . join(",",@escaped_IDs). "))";
    print STDERR "query=$query\n" if $debug;
    if ($link =~ /feature_id|patric_id/ and $opt->featurefields) { # get feature annotation
        my $featureFields = $opt->featurefields;
        $featureFields =~ s/ /_/g; # replace spaces with underscore
        $featureFields .= ",genome_id" if $opt->genomefields and $featureFields !~ /genome_id/;
        my $select = "select($featureFields,$link)"; 
        print STDERR "select=$select\n" if $debug;
        my %genome_to_links;
        my ($resp, $data) = $api->submit_query('genome_feature', "$query&$select&$limit");
        for my $record (@$data) {
            #print join("||", values %$record), "\n" if $debug;
            my $id = $record->{$link};
            print "feature $id\n" if $debug;
            for my $key (keys %$record) {
                next if $key eq $link;
                if ($featureFields =~ /$key/) {
                    print "     $key -> $record->{$key}\n" if $debug;
                    my $value = $record->{$key};
                    if (ref($value) eq 'ARRAY') {
                        $value = join(",", @$value); # flatten array to a string with all entries
                    }
                    $meta_column{$key}{$id} = $value;
                }
                if ($key eq 'genome_id') {
                    push(@{$genome_to_links{$record->{genome_id}}}, $id); # accumulate all links per genome_id
                }
            }
        }
        if (exists $meta_column{'genome_id'}) {
            my $genomeFields = $opt->genomefields();
            $genomeFields =~ s/ /_/g; # replace spaces with underscore
            if ($genomeFields) {
                my $unique_genomes = join(",", sort keys %genome_to_links);
                $query = "in(genome_id,($unique_genomes))";
                $select = "select($genomeFields,genome_id)";
                print STDERR "query db for genome fields:\n$query\n$select\n" if $debug;
                my ($resp, $data) = $api->submit_query('genome', "$query&$select&$limit");
                for my $record (@$data) {
                    #print join("||", keys %$record), "\n" if $debug;
                    my $genome_id = $record->{genome_id};
                    print "genome $genome_id\n" if $debug;
                    for my $key (keys %$record) {
                        if ($genomeFields =~ /$key/) {
                            print "  field $key: " if $debug;
                            for my $tree_id (@{$genome_to_links{$genome_id}}) {
                                my $value = $record->{$key};
                                if (ref($value) eq 'ARRAY') {
                                    $value = join(",", @$value); # flatten array to a string with all entries
                                }
                                $meta_column{$key}{$tree_id} = $value;
                                print " $tree_id => $record->{$key} " if $debug;
                            }
                            print "\n" if $debug;
                        }
                    }
                }
            }
        }
    }
    elsif ($link eq 'genome_id' and $opt->genomefields) { # get genome annotation
        my $fields = $opt->genomefields();
        $fields =~ s/ /_/g; # replace spaces with underscore
        my $select = "select($fields,genome_id)";
        if ($opt->verbose) {
            print STDERR "query = $query&$select&$limit\n\n";
        }
        my ($resp, $data) = $api->submit_query('genome', "$query&$select&$limit");
        print STDERR "data retrieved = $data\n" if $opt->verbose();
        for my $record (@$data) {
            print join("||", keys %$record)."\n" if $opt->verbose();
            my $id = $record->{genome_id};
            for my $key (keys %$record) {
                if ($fields =~ /$key/) {
                    my $value = $record->{$key};
                    if (ref($value) eq 'ARRAY') {
                        $value = join(",", @$value); # flatten array to a string with all entries
                    }
                    $meta_column{$key}{$id} = $value;
                }
            }
        }
    }
    for my $field (keys %meta_column) {
        $tree->add_tip_annotation($field, $meta_column{$field});
    }
}

if (0) {
    for my $column (sort keys %meta_column) {
    print "tree->add_phyloxml_tip_properties for $column\n" if $debug;
    my $provenance = "BVBRC";
    $provenance = $opt->provenance if $opt->provenance;
    $provenance = $meta_column{$column}{provenance} if exists $meta_column{$column}{provenance};
    $tree->add_tip_phyloxml_properties($meta_column{$column}, $column, $provenance);
}}

my $phyloxml_file = 'output.phyloxml';
if ($opt->output_name) {
    $phyloxml_file = $opt->ouptut_name;
    if ($phyloxml_file !~ /.phyloxml$/) {
        $phyloxml_file .= ".phyloxml";
    }
}
else {
    $phyloxml_file = $newickFile;
    $phyloxml_file =~ s/\.nwk$//;
    $phyloxml_file =~ s/\.tree$//;
    if ($opt->name_with_all_fields and $opt->databaselink) { # elaborate output file name with database fields added
        my $field_string = $opt->genomefields;
        if ($opt->databaselink ne 'genome_id') {
            $field_string .= "_" . $opt->featurefields;
        }
        $field_string =~ tr/,/_/;
        $phyloxml_file .= "_$field_string";
    }
    $phyloxml_file .= ".phyloxml";
}
open F, ">$phyloxml_file";
my $phyloxml_data = $tree->write_phyloXML();

if ($opt->remove_substring) {
    my $remove_substring = $opt->remove_substring;
    $phyloxml_data =~ s/$remove_substring//g;
}

print F $phyloxml_data;
close F;

if ($workspace_dir) {
    # copy phyoxml file to user's workspace
    my @command = ("p3-cp", "-m", "phyloxml=phyloxml");
    push(@command, " -f") if $opt->overwrite;
    push(@command, $phyloxml_file, "ws:" . $workspace_dir);
    print STDERR "commnd to copy back to workspace:\n@command\n" if $opt->verbose;
    my $rc = system(@command);
    die "Failure $rc running @command" unless $rc == 0;
    chdir($original_wd); # change back to the starting working directory
}
