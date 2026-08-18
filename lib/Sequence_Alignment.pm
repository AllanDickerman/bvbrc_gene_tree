package Sequence_Alignment;
use strict;
use warnings;
use List::Util ('any');
use P3DataAPI;
use IPC::Run 'run';

our $debug = 0;

sub new {
    my ($class, $input) = @_;
    if ($debug) {
        print(STDERR "in Sequence_Alignment::new\n");
    }
    my $self = {};
    bless $self, $class;
    $self->{_seqs} = {}; # allow multiple loci
    $self->{_ids} = [];
    $self->{_loci} = [];
    $self->{_default_locus} = '';
    $self->{_length} = {};
    $self->{_metadata} = {};

    if ($input) {
        print(STDERR " input file = $input\n") if $debug;
        $self->read_file($input)
    }
    print( STDERR "new Sequence_Alignment\n") if $debug;
    print(STDERR "loci = " . join(", ", @{$self->{_loci}}) . "\n") if $debug;
    return $self;
}

sub set_debug {my $onoff = shift; $debug = $onoff ? 1 : 0}
sub get_ntaxa { my $self = shift; return scalar(@{$self->{_ids}})}
sub get_length { 
    my ($self, $locus) = @_; 
    if (any {$locus eq $_} @{$self->{_loci}} ) {
        return($self->{_length}{$locus})
    }
    else {
        my $total_length = 0;
        for my $locus (@{$self->{_loci}}) {
            $total_length += $self->{_length}{$locus}
        }
        return $total_length;
    }
}
sub get_ids { 
    my ($self, $locus) = @_; 
    if ($locus) {
        return keys(%{$self->{_seqs}{$locus}})
    }
    return @{$self->{_ids}}
}

sub get_locus_ids {
    my $self = shift;
    return @{$self->{_loci}};
}

sub get_sequence_lengths {
    # returns hashref of hashrefs: {id}{locus} = length of sequence currently stored for that id/locus
    # reflects whatever state the object is in when called (raw or aligned)
    my $self = shift;
    my %lengths;
    for my $locus (@{$self->{_loci}}) {
        for my $id (keys %{$self->{_seqs}{$locus}}) {
            $lengths{$id}{$locus} = length($self->{_seqs}{$locus}{$id});
        }
    }
    return \%lengths;
}

sub add_seq { 
    my ($self, $id, $seq, $locus) = @_;
    if ($locus) {
        if (!$self->{_default_locus}) {
            $self->{_default_locus} = $locus;
            print  STDERR "setting default locus to '$locus'\n" if $debug;
        }    
    }
    else {
        if (!$self->{_default_locus}) {
            print STDERR "creating 'default' locus name\n" if $debug;
            $self->{_default_locus} = 'default';
        }
        $locus = $self->{_default_locus};
    }
    print("add_seq $id " . substr($seq, 0, 10) . " locus= '$locus') \n") if $debug;
    if (! exists $self->{_seqs}{$locus}) {
        $self->{_seqs}{$locus} = {};
        $self->{_length}{$locus} = 0;
        push @{$self->{_loci}}, $locus;
        my @sorted_loci = sort(@{$self->{_loci}});
        $self->{_loci} = \@sorted_loci;
    }
    if ($self->{_in_alignment}) {
        # if _in_alignment we are reading aligned version of sequence, just replace existing sequence
        print STDERR "reading aligned version of $id at locus $locus length = ", length($seq), "\n" if $debug;
    }
    else {
        print STDERR "adding id $id to locus $locus\n" if $debug;
        if (exists $self->{_seqs}{$locus}{$id}) { # make identifier unique in case of duplicate
            warn("duplicate occurence of id $id at locus $locus");
            my $temp = $id;
            my $suffix = 1;
            while (exists $self->{_seqs}{$locus}{$temp}) {
                $suffix++;
                $temp = "${id}_$suffix";
            }
            print STDERR "sequence id $id exists\nincrementing to $temp\n" if $debug;
            $id = $temp;
        }
    }
    $seq =~ tr/ //d;
    $self->{_seqs}{$locus}{$id} = $seq;
    if (length($seq) > $self->{_length}{$locus}) {
        $self->{_length}{$locus} = length($seq)
    }
    # add id to _ids if it is not already there (keep it unique)
    unless (grep(/^$id$/, @{$self->{_ids}})) {
        push(@{$self->{_ids}}, $id);
    }
    return $id; # might be different from input id
}

sub set_metadata {
    # allow stroring abritrary data per id
    my ($self, $id, $category, $value) = @_;
    if (!exists $self->{_metadata}{$category}) {
        $self->{_metadata}{$category} = {};
    }
    $self->{_metadata}{$category}{$id} = $value;
}

sub get_metadata {
    my ($self, $id, $category) = @_;
    if (exists $self->{_metadata}{$category}) {
        return $self->{_metadata}{$category}{$id};
    }
    return undef;
}

sub detect_file_format {
    my $class = shift;
    my $fh = shift;
    print STDERR "in detect_format, class=$class, fh=$fh\n" if $debug;

    $_ = readline $fh;
    seek $fh, 0, 0; # reset file to beginning

    print STDERR "first line of file is :\n", $_, "\n" if $debug;
    warn "cannot read first line of file" unless $_;
    my $format = 'unknown';
    $format = 'clustal' if (/^CLUSTAL/ || /^MUSCLE/);
    $format = 'fasta' if (/^>/);
    $format = 'phylip' if (/^(\d+)\s+(\d+)\s*$/);
    $format = 'nexus' if (/\#NEXUS/);
    print STDERR "input format detected as $format\n" if $debug;
    return $format
}

sub read_file {
    my ($self, $fh, $format, $locus) = @_;
    my @retval; #return seq ids
    my $filename = $fh;
    if ( ! ref($fh) ) {
        $fh = undef;
        open $fh, $filename;
        if (! $locus) {
            $locus = $filename;
            $locus =~ s/.*\///; #delete up to last forward slash
            $locus =~ s/\..*//; #delete after last dot
        }
    } 
    if (! $locus) {
        $locus = $self->{_default_locus};
    }
    if (! $format) {
        $format = $self->detect_file_format($fh);
        print("format detected as $format") if $debug;
    }
    print STDERR "read_file $fh, $format, $locus \n" if $debug;
    if ($format eq 'unknown') {
        return undef;
    }
    if ($format eq 'clustal') {
        my $found = 0;
        while (<$fh>) {
            if (/^CLUSTAL/ || /^MUSCLE/) {
                $found = 1;
                last;
            }
        }
        die "Format seems to be wrong, not Clustal.\n" if (!$found); 
        while (<$fh>) {
            if (/^(\S+)\s+(\S+)/) {
                my ($id, $seq) = ($1, $2);
                $id = $self->add_seq($id, $seq, $locus);
                push @retval, $id;
            }
        }
    }
    elsif ($format eq 'phylip') {
        $_ = <$fh>;
        die "Format does not seem to be phylip\n" unless (/^\s*(\d+)\s+(\d+)\s*$/);
        my $ntaxa = $1;
        my $nchar = $2;
        my %seqhash = {};
        my @ids;
        for my $i (1..$ntaxa) {
            $_ = <$fh>;
            /(\S+)\s+(\S.*\S)/ or die $_;
            my $id = $1;
            my $seq = $2;
            push(@ids, $id);
            $seqhash{$id} = $seq;
        }
        # now if there are more lines, read in same order as first set, but without identifiers
        my $index = 0;
        while (<$fh>) {
            my $seq = $_;
            $seq =~ s/\s//g;
            if ($seq) {
                my $id = $ids[$index % $ntaxa];
                $seqhash{$id} .= $seq;
                $index += 1
            }
        }
        foreach my $id (@ids) { # phylip uses '.' as insert (unknown) character
            $seqhash{$id} =~ s/\./\-/g; # replace dot as gap char with '-'
            $id = $self->add_seq($id, $seqhash{$id}, $locus);
            push @retval, $id;
        }
    }
    elsif ($format eq 'fasta') {
        #while (    my($id, $def, $seq) = gjoseqlib::read_next_fasta(\*$fh))
        my ($id, $seq);
        while (<$fh>) {
            if (/^>(\S+)/) {
                if ($seq) {
                    $id = $self->add_seq($id, $seq, $locus);
                    push @retval, $id;
                    $seq = '';
                }
                $id = $1;
            }
            elsif (/(.*\S)/) {
                $seq .= $1;
            }
        }
        if ($seq) {       
            $id = $self->add_seq($id, $seq, $locus);
            push @retval, $id;
        }
    }
    elsif ($format eq 'nexus')
    {
        $_ = <$fh>;
        die "Format does not seem to be NEXUS" unless (/\#NEXUS/);
        my $matrix;
        while (<$fh>) {
            chomp;
            s/\[[^\]]*\]//g; # strip out comments in square brackets
            if (/^matrix/i){
                $matrix = 1;
                next;
            }
            if ($matrix) {	    
                if (/^(\S+)\s+(\S+)/) {
                    my ($id, $seq) = ($1, $2);
                    $id = $self->add_seq($id, $seq, $locus);
                    push @retval, $id;
                }
                last if (/;/);
            }
        }
    }
    print STDERR "Done reading $filename for locus $locus, len=$self->{_length}{$locus}, num=", scalar(@{$self->{_loci}}), "\n" if $debug;
    print STDERR "Number of ids = ", scalar @{$self->{_ids}}, "\n" if $debug;
    my @loci = $self->get_locus_ids();
    print STDERR "number of loci is ", scalar(@loci), " ", join(',', @loci), "\n" if $debug;
    return @retval;
}

sub write_to_file {
    my ($self, $fh, $format, $locus) = @_;
    my $out;
    if ( ref($fh) ) {
        $out = $fh;
    }
    else {
        print STDERR "in write_to_file, not a file handle, open file $fh\n" if $debug;
        open($out, ">", $fh);
    } 
    if (! $format) {
        $format = 'fasta';
    }
    my @ids = $self->get_ids();
    print STDERR "number of ids is ", scalar(@ids), "\n" if $debug;
    my %seqs = ();
    my @loci;
    if (any { $_ eq $locus} @{$self->{_loci}}) {
        @loci = ($locus);
        print "limit to locus $locus\n" if $debug;
    }
    else {
        @loci = $self->get_locus_ids();
        print STDERR "concatenate all ", scalar(@loci), " loci: ", join(',', @loci), "\n" if $debug;
    }
    for my $loc (@loci) {
        print STDERR "ids at locus $loc: ", join(',', keys(%{$self->{_seqs}{$loc}})), "\n" if $debug;
        for my $id (@ids) {
            if (exists($self->{_seqs}{$loc}{$id}) and $self->{_seqs}{$loc}{$id}) {
                $seqs{$id} .= $self->{_seqs}{$loc}{$id};
                print STDERR "$id $loc ", length($self->{_seqs}{$loc}{$id}), "   ",substr($self->{_seqs}{$loc}{$id},0,23), "\n" if $debug;
            }
            else {
                print(STDERR "id $id does not exist at locus $loc, use " . $self->{_length}{$loc} . " gaps\n");
                $seqs{$id} .= '-'x$self->get_length($loc);
            }
        }
    }
    my $max_id_length = 0;
    if ($format eq 'phylip') {
        print $out $self->get_ntaxa(), "  ", $self->get_length($locus), "\n";
        for my $id (@ids) {
            if (length($id) > $max_id_length) {
                $max_id_length = length($id);
            }
        }
    }
    for my $id (@ids) {
        print STDERR "writing $id, length = ",length($seqs{$id}), "\n" if $debug; 

        if (($format eq 'fasta') or ($format eq 'raxml')) {
            if ($self->{_in_alignment}) { # writing unaligned seqs for purpose of aligning
                $seqs{$id} =~ tr/-//d;
            }
            my $output_id = $id;
            if ($format eq 'raxml') {
                my $changed = $output_id =~ tr/:()[]/_____/; #replace with underscores
                if ($changed) {
                    print STDERR "for raxml: original=$id, changed=$output_id\n" if $debug;
                    if (!exists( $self->{_raxml_to_original_id})) {
                        $self->{_raxml_to_original_id} = {};
                    }
                    $self->{_raxml_to_original_id}{$output_id} = $id;
                }
            }
            print $out ">$output_id\n$seqs{$id}\n";
        }
        elsif ($format eq 'phylip') {
            #print $out $id, " "x($max_id_length-length($id)+4), $seqs{$id}, "\n";
            printf($out "%-${max_id_length}s  %s\n", $id, $seqs{$id});
        }
    }
    if ($out ne $fh) {
        print "closing $fh\n" if $debug;
        close $out  #because we opened it
    } 
}

sub write_fasta {
  # write out in fasta format
    my ($self, $out, $locus) = @_;
    $self->write_to_file($out, 'fasta', $locus);
}

sub write_phylip {
  # write out in phylip format
    my ($self, $out, $locus) = @_;
    $self->write_to_file($out, 'phylip', $locus);
}

sub write_fasta_for_raxml { # edit out illegal characters in sequence ids
    my ($self, $out, $locus) = @_;
    $self->write_to_file($out, 'raxml', $locus);
}

sub restore_original_ids_in_raxml_tree {
    my ($self, $newick) = @_;
    return $newick unless exists $self->{_raxml_to_original_id};
    for my $raxml_id (keys %{$self->{_raxml_to_original_id}}) {
        $newick =~ s/$raxml_id/$self->{_raxml_to_original_id}{$raxml_id}/;
    }
    return $newick
}

sub align {
    my ($self, $locus) = @_;
    my @command_list = ();
    my $details = "";
    if ($locus) {
        print STDERR "align($locus)\n" if $debug;
        my $unaligned = "${locus}_unaligned.fa";
        $self->write_to_file($unaligned, "fasta", $locus);
        my $parallel = $ENV{P3_ALLOCATED_CPU};

        my $cmd = ["mafft", "--auto"];
        if ($parallel) {
            push @$cmd, "--thread", $parallel;
        }
        push @$cmd, $unaligned;
        my $comment = join(" ", @$cmd);
        print STDERR $comment, "\n" if $debug;
        my $aligned_fasta = "${locus}_aligned.fa";
        print STDERR  "Align with mafft\n" if $debug;
        my $stderr = "alignment_stderr.txt";
        unlink($stderr);
        my $retval = run($cmd, ">", $aligned_fasta, "2>", $stderr);
        $self->{_in_alignment} = 1;
        $self->read_file($aligned_fasta, 'fasta', $locus);
        delete($self->{_in_alignment});
        open F, $stderr;
        $details = join(" ", @$cmd)."\n";
        $details .= do {
            local $/; # Undefines the line separator locally
            <F>;    # Reads the whole file as a single string
        };
        close F;
        push @command_list, join(" ", @$cmd)."\n";
        print STDERR "command_list = ", join("\n\ncmd=", @command_list), "\n";
    }
    else {
        print STDERR "align()\n" if $debug;
        for my $locus (@{$self->{_loci}}) {
            my ($cmd_aref, $detail) = $self->align($locus);
            push @command_list, $cmd_aref->[0];
            $details .= $detail;
        }
    }
    return \@command_list, $details; # return ref to list of commands and a block of text with details (stderr)
}

sub calc_column_gap_count {
    my ($self, $locus) = @_;
    print STDERR "calc_column_gap_count($locus)\n";
    my @gap_count;
    $#gap_count = $self->{_length}{$locus}-1;
    for my $id (keys %{$self->{_seqs}{$locus}}) {
        my @str_as_array = split('', $self->{_seqs}{$locus}{$id});
        for my $gap_pos (0 .. $#str_as_array) {
            $gap_count[$gap_pos] += $str_as_array[$gap_pos] eq '-';
        }
    }
    if ($debug & 0) {
        my $vis = '';
        print "Length of \@gap_count = ", scalar(@gap_count), ", vs self->length = $self->{_length}{$locus}\n";
        for my $i (0 .. $self->{_length}{$locus}-1 ) {
            my $prop10 = int(10 * $gap_count[$i] / $self->get_ntaxa());
            $prop10 = 9 if $prop10 > 9;
            $vis .= $prop10;
        }
        print "$vis\n";
    }
    return \@gap_count;
}

sub end_trim {
    # trim gappy ends inward to a minimum occupancy threshold (proportion of non-gap chars)
    my ($self, $threshold, $locus) = @_;
    ($threshold <= 1.0 and $threshold > 0) or die "threshold must be between 0 and 1";
    my @steps;
    my $details = "";
    if ($locus) {
        my $comment = sprintf("Evaluate $locus alignment, trim ends over %.0f percent gaps.\n", $threshold*100);
        #push @steps, $comment;
        $details .= $comment;
        print STDERR $comment;
        print STDERR "end_trim($threshold, $locus)\n";
        my $gap_count = $self->calc_column_gap_count($locus);
        my $max_gaps = (1.0 - $threshold) * scalar(keys %{$self->{_seqs}{$locus}});
        print STDERR "max_gaps = $max_gaps\n" if $debug;
        my $total_gaps = 0;
        for my $i (0 .. $self->{_length}{$locus}-1 ) {
            $total_gaps += $gap_count->[$i];
        }
        my $prop_gaps = $total_gaps / ($self->{_length}{$locus} * scalar keys %{$self->{_seqs}{$locus}});
        $comment = sprintf("Before trimming, length = $self->{_length}{$locus}, proportion gaps = %.4f\n", $prop_gaps);
        $details .= $comment;
        print STDERR $comment;

        my $left_trim = 0;
        $left_trim++ while (($gap_count->[$left_trim] > $max_gaps) and $left_trim < $self->{_length}{$locus}-1);
        my $right_trim = 0;
        $right_trim++ while (($gap_count->[$self->{_length}{$locus} - 1 - $right_trim] > $max_gaps) and $right_trim < $self->{_length}{$locus}-1);
        $right_trim++ while (($right_trim > $left_trim) and ($gap_count->[$right_trim-1] > $max_gaps));
        print(STDERR "trim $left_trim positions at left\n");
        print(STDERR "trim $right_trim positions at right\n");
        if ($left_trim or $right_trim) {
            my $len = $self->{_length}{$locus} - $right_trim - $left_trim;
            $comment = "Trim $left_trim on left, trim $right_trim on right, new length = $len.\n";
            push @steps, $comment;
            $details .= $comment;
            print STDERR $comment; 
            for my $id (keys %{$self->{_seqs}{$locus}}) {
                $self->{_seqs}{$locus}{$id} = substr($self->{_seqs}{$locus}{$id}, $left_trim, $len);
            }
            $self->{_length}{$locus} = $len;
            $total_gaps = 0;
            for my $i ($left_trim .. scalar(@$gap_count)-$right_trim-1 ) {
                $total_gaps += $gap_count->[$i];
            }
            $prop_gaps = $total_gaps / ($self->{_length}{$locus} * scalar keys %{$self->{_seqs}{$locus}});
            $comment = sprintf("After trimming, length = $self->{_length}{$locus}, proportion gaps = %.4f\n", $prop_gaps);
            $details .= $comment;
            print(STDERR $comment);
        }
        else { 
            $comment = "No trimming needed for $locus\n";
            push @steps, $comment;
            $details .= $comment;
        } 
        print STDERR "At end of end_trim($threshold, $locus), step = \n" . join(" ", @steps);
    }
    else { # locus not passed
        print STDERR "In end_trim($threshold)\n" if $debug;
        my $details;
        for my $locus (@{$self->{_loci}}) {
            my ($locus_steps, $locus_detail) = $self->end_trim($threshold, $locus);
            push @steps, @{$locus_steps};
            $details .= $locus_detail;
            print STDERR "";
        }
    }
    return \@steps, $details; # return ref to list of commands and a block of text with details (stderr)
}

sub calc_row_gap_count {
    my ($self, $locus) = @_;
    print STDERR "In calc_row_gap_count($locus)\n" if $debug;
    my %gap_count;
    for my $id (keys %{$self->{_seqs}{$locus}}) {
        $gap_count{$id} = $self->{_seqs}{$locus}{$id} =~ tr/-/-/;
    }
    return \%gap_count;
}

sub get_sequence_occupancy {
    # returns hashref: {id}{locus} = percentage of non-gap characters in that id's aligned sequence at that locus
    # ids lacking a sequence at a given locus simply have no entry there (caller can detect via exists)
    my $self = shift;
    my %occupancy;
    for my $locus (@{$self->{_loci}}) {
        my $length = $self->{_length}{$locus};
        my $gaps_per_seq = $self->calc_row_gap_count($locus);
        for my $id (keys %$gaps_per_seq) {
            $occupancy{$id}{$locus} = $length ? 100 * (1 - $gaps_per_seq->{$id} / $length) : 0;
        }
    }
    return \%occupancy;
}

sub delete_gappy_seqs {
    # remove gappy sequences below minimum occupancy threshold (proportion of non-gap chars)
    my ($self, $threshold, $locus) = @_;
    ($threshold <= 1.0 and $threshold > 0) or die "threshold must be between 0 and 1";
    my @retval;
    if ($locus) {
        print STDERR "In delete_gappy_seqs($threshold, $locus)\n" if $debug;
        my $max_gaps = (1.0 - $threshold)*$self->get_length($locus);
        for my $id (keys %{$self->{_seqs}{$locus}}) {
            my $num_gaps = $self->{_seqs}{$locus}{$id} =~ tr/-/-/;
            if ($num_gaps > $max_gaps) {
                delete($self->{_seqs}{$locus}{$id});
                my $command = "sequence $id";
                $command .= " at segment $locus" unless $locus eq 'default';
                $command .= sprintf(" deleted due to $num_gaps gaps (%.2f%%), exceeding threshold of %.2f%%\n", $num_gaps/$self->get_length($locus), $threshold);
                print STDERR $command;
                push @retval, $command;
            }
        }
        if (scalar @retval) {
            print STDERR "recalc locus length\n" if $debug;
            $self->{_length}{$locus} = 0;
            for my $id (keys %{$self->{_seqs}{$locus}}) {
                if (length($self->{_seqs}{$locus}{$id}) > $self->{_length}{$locus}) {
                    $self->{_length}{$locus} =  length($self->{_seqs}{$locus}{$id});
                }
            }
        }
    }
    else {
        print STDERR "In delete_gappy_seqs($threshold)\n" if $debug;
        for my $locus (@{$self->{_loci}}) {
            my $deleted_seqs = $self->delete_gappy_seqs($threshold, $locus);
            if ($deleted_seqs and scalar(@$deleted_seqs)) {
                push @retval, @$deleted_seqs;
            }
        }
        unless (scalar @retval) {
           @retval = ("No sequences gappy enough to remove");
        }
    }
    return \@retval;
}

sub write_stats {
    my $self = shift;
    print STDERR "write_stats()\n" if $debug;
    my %prop_gaps;
    for my $locus (@{$self->{_loci}}) {
        my $gaps_per_seq = $self->calc_row_gap_count($locus);
        my $total_gaps = 0;
        for my $id (keys %$gaps_per_seq) {
            $total_gaps += $gaps_per_seq->{$id};
        }
        $prop_gaps{$locus} = $total_gaps / ((scalar keys %$gaps_per_seq) * $self->{_length}{$locus});
    }

    my $retval = "";
    $retval .= "Alignment Statistics   " . join("\t", @{$self->{_loci}}) . "\n";
    $retval .= "Number of sequences: ";
    for my $locus (@{$self->{_loci}}) {
        $retval .= "\t" . scalar(keys %{$self->{_seqs}{$locus}});
    }
    $retval .= "\n";

    $retval .= "Alignment length:    ";
    for my $locus (@{$self->{_loci}}) {
        $retval .= "\t" . $self->get_length($locus);
    }
    $retval .= "\n";

    $retval .= "Proportion gaps:     ";
    for my $locus (@{$self->{_loci}}) {
        $retval .= "\t" . sprintf("%.4f", $prop_gaps{$locus});
    }
    $retval .= "\n";


    return $retval;
}


return 1
