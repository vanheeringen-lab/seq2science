#!/usr/bin/perl
use strict;
use warnings;

my $num_args = $#ARGV + 1;
if ($num_args != 2) {
    print "Usage: atacshift.pl input.sam genome.fa.sizes\n";
    exit;
}

open my $in_sizes, "<:encoding(utf8)", $ARGV[1];

my %sizes;
while (my $line = <$in_sizes>) {
  chomp $line;

  my @split_line = split '\t', $line;
  $sizes{$split_line[0]} = $split_line[1];
}

open my $in_sam, "<:encoding(utf8)", $ARGV[0];

while (my $line = <$in_sam>) {
  chomp $line;

  # ignore header lines
  if ($line =~ (/^(\@)/)) {
    print "$line\n";
    next;
  }

  # split line into array
  my @read = split '\t', $line;

  # only keep (necessary) first 11 columns
  @read = @read[0..10];
  
  if (not ($read[1] & 0x04)) {  # read must be mapped to shift
    # my ($read_len, $_mapping) = split(/[A-Z]/, $read[5], 2);
    #length($read[10]);
    
    # empty QUAL and SEQ
    $read[9] = "*";
    $read[10] = "*";

    my @cigar = ( $read[5] =~ /\d+M/g );
    my $read_len;
    if ($read[1] & 0x10) {  # reverse strand
      if (($read[1] & 0x01) && not ($read[1] & 0x08)) {  # read paired and mate not unmapped
          $read[7] = $read[7] + 4;  # mate shifts 4
          $read[8] = $read[8] + 9;  # template length decreases 4 + 5
      }
      # $read_len = (split /[A-Z]/, $read[5])[-1];
      $read_len = substr $cigar[-1], 0, -1;
      $read_len -= 5;
    }
    else {  # forward strand
      if (($read[1] & 0x01) && not ($read[1] & 0x08)) {  # read paired and mate not unmapped
          $read[8] = $read[8] - 9;  # template length decreases 4 + 5
      }
      $read[3] = $read[3] + 4;  # shift the position +4
      #($read_len, $_mapping) = split(/[A-Z]/, $read[5], 2);
      $read_len = substr $cigar[0], 0, -1;
      $read_len -= 4;
    }
    $read[5] = "${read_len}M";
  }
  if (($read[3] < 0) || ($read[3] >= $sizes{$read[2]})) {
    print STDERR "MAMA MIA $read[3]";
    next;
  }
  print join "\t", @read, "\n";
}
