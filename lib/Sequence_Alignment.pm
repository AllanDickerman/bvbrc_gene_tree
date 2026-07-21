package Sequence_Alignment;
use strict;
use warnings;
use List::Util qw(max);
#
# Don't import all; we are clashing with write_fasta
#use gjoseqlib qw();

our $debug = 0;

sub new {
    my ($class, $input) = @_;
    if ($debug) {
        print(STDERR "in Sequence_Alignment::new\n");
        print(STDERR " class = $class\n");
        print(STDERR " input = $input\n");
    }
    my $self = {};
    bless $self, $class;
    $self->{_seqs} = {}; # allow multiple loci
    #$self->{_annot} = {};
    $self->{_ids} = [];
    $self->{_default_locus} = '';
    #    $self->{_is_aligned} = 0;
    $self->{_length} = {};

    if ($input) {
        $self->read_file($input)
    }
    print( STDERR "new Sequence_Alignment\n") if $debug;
    return $self;
}

sub set_debug {my $onoff = shift; $debug = $onoff ? 1 : 0}
#sub set_is_aligned {my $self = shift; my $onoff = shift; $self->{_is_aligned} = $onoff ? 1 : 0}
sub get_ntaxa { my $self = shift; return scalar(@{$self->{_ids}})}
sub get_length { 
    my ($self, $locus) = @_; 
    if ($locus) {
        return($self->{_length}{$locus})
    }
    else {
        my $total_length = 0;
        for my $locus (keys %{$self->{_length}}) {
            $total_length += $self->{_length}{$locus}
        }
        return $total_length;
    }
}
sub get_ids { 
    my ($self, $locus) = @_; 
    if ($locus) {
        return keys(%{$self->{_seqs}})
    }
    return $self->{_ids}
}

sub add_seq { 
    my ($self, $id, $seq, $locus) = @_;
    if ($locus) {
        if (!$self->{_default_locus}) {
            $self->{_default_locus} = $locus;
        }    
    }
    else {
        if (!$self->{_default_locus}) {
            $self->{_default_locus} = 'default';
        }
        $locus = $self->{_default_locus};
    }

    if (! exists $self->{_seqs}{$locus}) {
        $self->{_seqs}{$locus} = {};
        $self->{_length}{$locus} = 0;
    }
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
    $self->{_seqs}{$locus}{$id} = $seq;
    if (length($seq) > $self->{_length}{$locus}) {
        $self->{_length}{$locus} = length($seq)
    }
    # add id to _ids if it is not already there (keep it unique)
    unless (exists($self->{_ids}[$id])) {
        push(@{$self->{_ids}}, $id);
    }
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
    if ( ! ref($fh) ) {
        print STDERR "in read_file, not a file handle, open file $fh\n" if $debug;
        my $temp = undef;
        open $temp, $fh;
        $fh = $temp;
    } 
    if (! $locus) {
        $locus = $self->{_default_locus};
    }
    if (! $format) {
        $format = $self->detect_file_format($fh)
    }
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
                $self->add_seq($id, $seq, $locus)
            }
        }
    }
    elsif ($format eq 'phylip') {
        $_ = <$fh>;
        die "Format does not seem to be phylip\n" if (!/^\s*(\d+)\s+(\d+)\s*$/);
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
            $self->add_seq($id, $seqhash{$id}, $locus);
        }
    }
    elsif ($format eq 'fasta') {
        #while (    my($id, $def, $seq) = gjoseqlib::read_next_fasta(\*$fh))
        #
        my ($id, $seq);
        while (<$fh>) {
            if (/^>(\S)+/) {
                if ($seq) {
                    $seq =~ tr/ //d;
                    $self->add_seq($id, $seq, $locus);
                }
                /^>(\S)+/;
                $id = $1;
            }
            elsif (/(.*\S)/) {
                $seq .= $1;
            }
        }
        if ($seq) {       
            $seq =~ tr/ //d;
            $self->add_seq($id, $seq, $locus)
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
                    $self->add_seq($id, $seq, $locus)
                }
                last if (/;/);
            }
        }
    }
}

sub write_to_file {
    my ($self, $fh, $format, $locus) = @_;
    if (!$locus) {
        $locus = $self->{_default_locus};
    }
    my $out = $fh;
    if ( ! ref($fh) ) {
        print STDERR "in write_to_file, not a file handle, open file $fh\n" if $debug;
        open($out, ">", $fh);
    } 
    if (! $format) {
        $format = 'fasta';
    }
    my @ids = $self->get_ids($locus);
    my %seqs;
    my @loci = [$locus];
    if ($locus eq 'all') {
        @loci = $self->get_locus_ids();
    }
    for my $loc (@loci) {
        for my $id (@ids) {
            if (exists($self->{_seqs}{$loc}{$id})) {
                $seqs{$id} .= $self->{_seqs}{$loc}{$id}
            }
            else {
                $seqs{$id} .= '-'*$self->get_length($loc)
            }
        }
    }
    my $max_id_length = 0;
    if ($format eq 'phylip') {
        print $out $self->get_ntaxa(), "  ", $self->get_length(), "\n";
        for my $id (@ids) {
            if (length($id) > $max_id_length) {
                $max_id_length = length($id);
            }
        }
    }
    if ($format eq 'raxml') {
        #my $raxml_illegal_chars = ":()[]";
        for my $id (@{$self->{_ids}}) {
            if ($id =~ tr/:()[]/:()[]/) { # counts but doesn't change
                $self->{_raxml_to_orignal_id} = {};
                last
            }
        }
    }
    for my $id (@ids) {
        if ($format eq 'fasta') {
            print($out, ">$id\n$seqs{$id}\n");
        }
        elsif ($format eq 'phylip') {
            print($out, $id, " "*($max_id_length-length($id)+4), $seqs{$id}, "\n");
            #printf($out "%-${max_id_length}s  %s\n", $id, $seq);
        }
        #elsif ($format eq 'nexus') {
        #    print($out, $id, " "*(max_id_length-length($id)+2), $seq{$id}, "\n");
        #}
        elsif ($format eq 'raxml') {
            my $orig = $id;
            my $seq = $seqs{$id};
            my $changed = $id =~ tr/:()[]/_____/; #replace with underscores
            if ($changed) {
                print STDERR "in write_fasta_for_raxml: original=$orig, changed=$id\n" if $debug;
                if (!exists( $self->{_raxml_to_orignal_id})) {
                    $self->{_raxml_to_orignal_id} = {};
                }
                $self->{_raxml_to_original_id}{$id} = $orig;
            }
            print $out ">$id\n$seq\n";
        }
    }
    if ($out ne $fh) {
        print "closing $fh\n" if $self->{_debug};
        close $out  #because we opened it
    } 
}

sub unalign() {
    my $self = shift;
    foreach my $locus (keys %{$self->{_seqs}}) {
        foreach my $id (keys %{$self->{_seqs}{$locus}}) {
                $self->{_seqs}{$locus}{$id} =~ tr/-//d;
        }
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
    $self->write_to_file($out, 'fasta', $locus);
}

sub write_fasta_for_raxml { # edit out illegal characters in sequence ids
    my ($self, $out, $locus) = @_;
    $self->write_to_file($out, 'raxml', $locus);
}

sub restore_original_ids_in_raxml_tree {
    my ($self, $newick) = @_;
    return $newick unless exists $self->{_raxml_to_orignal_id};
    for my $raxml_id (keys %{$self->{_raxml_to_orignal_id}}) {
        $newick =~ s/$raxml_id/$self->{_raxml_to_orignal_id}{$raxml_id}/;
    }
    return $newick
}

sub calc_column_gap_count {
    my ($self, $locus) = shift;
    #print STDERR "In calc_column_gap_count\n";
    my @gap_count;
    $#gap_count = $self->{_length}{$locus}-1;
    for my $id (keys %{$self->{_seqs}{$locus}}) {
        my @str_as_array = split('', $self->{_seqs}{$locus}{$id});
        for my $gap_pos (0 .. $#str_as_array) {
            $gap_count[$gap_pos] += $str_as_array[$gap_pos] eq '-';
        }
    }
    return \@gap_count;
}

sub write_stats {
    my $self = shift;
    my $gaps_per_seq = $self->calc_row_gap_count();
    my $worst_seq = undef;
    my $worst_seq_gaps = 0;
    my $avg_gaps_per_seq = 0;
    my $total_gaps = 0;
    for my $id (keys %$gaps_per_seq) {
        $avg_gaps_per_seq += $gaps_per_seq->{$id};
        $total_gaps += $gaps_per_seq->{$id};
        if ($gaps_per_seq->{$id} > $worst_seq_gaps) {
            $worst_seq_gaps = $gaps_per_seq->{$id};
            $worst_seq = $id
        }
    }
    $avg_gaps_per_seq /= $self->get_ntaxa();

    my $retval = "";
    $retval .= "Alignment Statistics\n";
    $retval .= "\tNumber of sequences    = ". $self->get_ntaxa(). "\n";
    $retval .= "\tAlignment length       = ". $self->get_length(). "\n";
    $retval .= "\tProportion gaps        = ". sprintf("%.4f", $avg_gaps_per_seq/$self->get_length()). "\n";
    #$retval .= "\tAverage column entropy = ". sprintf("%.3f", $avg_entropy). "\n";
    return $retval;
}

sub end_trim {
    # trim gappy ends inward to a minimum occupancy threshold (proportion of non-gap chars)
    my ($self, $threshold, $locus) = @_;
    print STDERR "In end_trim($threshold)\n" if $debug;
    ($threshold <= 1.0 and $threshold > 0) or die "threshold must be between 0 and 1";
    if ($locus) {
        my $gap_count = $self->calc_column_gap_count($locus);
        my $max_gaps = (1.0 - $threshold) * scalar(keys %{$self->{_seqs}{$locus}});
        if ($debug & 0) {
            my $vis = '';
            my $vis2 = '';
            print "Length of \@gap_count = ", scalar(@$gap_count), ", vs self->length = $self->{_length}\n";
            for my $i (0 .. $self->{_length}-1) {
                my $prop10 = int(10 * $gap_count->[$i] / $self->get_ntaxa());
                $prop10 = 9 if $prop10 > 9;
                $vis .= $prop10;
                $vis2 .= $i % 10;
            }
            print "$vis\n$vis2\n";
        }

        my $start = 0;
        $start++ while ($gap_count->[$start] > $max_gaps and $start < $self->{_length}-1);
        my $end = $self->{_length}-1;
        $end-- while ($end and $gap_count->[$end] > $max_gaps);
        my $num_end_columns_trimmed = $self->{_length} - $end - 1;
        my $len = $end - $start + 1;
        print STDERR "Trim up to $start and after $end\n" if $debug;
        #print substr($vis, $start, $len), "\n";
        #print substr($vis2, $start, $len), "\n";
        for my $id (keys %{$self->{_seqs}{$locus}}) {
                $self->{_seqs}{$locus}{$id} = substr($self->{_seqs}{$locus}{$id}, $start, $len);
                #die "end trimming made sequence $id too short: ", length($self->{_seqs}->{$id}), " vs $len" if length($self->{_seqs}->{$id}) != $len;
        }
        $self->{_length} = $len;
        return ($start, $num_end_columns_trimmed);
    }
    else { # locus not passed
        for my $locus (keys %{$self->{_seqs}}) {
            $self->end_trim($threshold, $locus)
        }
    }
}

sub calc_row_gap_count
    {
        my $self = shift;
        print STDERR "In calc_row_gap_count\n" if $debug;
        my %gap_count;
        for my $id (@{$self->{_ids}}) {
            $gap_count{$id} = $self->{_seqs}->{$id} =~ tr/-/-/;
        }
        return \%gap_count;
    }

    sub delete_gappy_seqs {
        # remove gappy sequences below minimum occupancy threshold (proportion of non-gap chars)
        my $self = shift;
        my $threshold = shift;
        print STDERR "In delete_gappy_seqs($threshold)\n" if $debug;
        ($threshold <= 1.0 and $threshold > 0) or die "threshold must be between 0 and 1";
        my $gap_count = $self->calc_row_gap_count();
    my $max_gaps = (1.0 - $threshold)*$self->{_length};
    my $index = 0;
    my @retval = ();
    while ($index <= $#{$self->{_ids}}) {
        my $id = $self->{_ids}->[$index];
        if ($id ne $self->{_ids}->[$index]) {
            print "item at $index is not $id: ", join(", ", @{$self->{_ids}}[$index - 1, $index, $index + 1]), "\n";
        }
        if ($gap_count->{$id} > $max_gaps) {
            # remove id from list and seq from hash
            print STDERR "seq $id has $gap_count->{$id} gaps, deleting index $index.\n" if $debug;
            my $hashlen_before = scalar keys %{$self->{_seqs}};
            my $arraylen_before = scalar @{$self->{_ids}};
            my @abefore= @{$self->{_ids}}[$index-1,$index, $index+1];
            delete($self->{_seqs}->{$id});
            my $removed = splice @{$self->{_ids}}, $index, 1;
            my $hashlen_after = scalar keys %{$self->{_seqs}};
            my $arraylen_after = scalar @{$self->{_ids}};
            if ($self->{_ids}->[$index] eq $id) {
                print "list before = ",join(", ", @abefore), "\n"; 
                print "list after  = ", join(", ", @{$self->{_ids}}[$index - 1, $index, $index + 1]), "\n";
                print "Removed = $removed\n";
                die "found $id on list after deleting it, index=$index hashlen_b=$hashlen_before, hashlen_a=$hashlen_after; arraylen_b=$arraylen_before, arraylen_a=$arraylen_after"; 
            }
            push @retval, $id;
        }
        else {
            $index++;
        }
    }
    die "array vs hash mismatch after deteting gappy seqs" if (scalar @{$self->{_ids}} != scalar keys %{$self->{_seqs}});
    return \@retval
}    


return 1
