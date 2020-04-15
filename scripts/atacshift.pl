#!/usr/bin/perl
#
# atacshift: Correct for tn5 bias by shifting reads on the positive strand +4 and reads 
#            on the negative strand -5. The CIGAR string is the same length as the read
#            but all nucleotides now map (M). Reads that fall (partly) outside of the 
#            genome after shifting are removed.
#  
# author: Maarten van der Sande

use strict;
use warnings;
use List::Util qw(sum);


my %sizes;
my $num_args = $#ARGV + 1;
if ($num_args != 2) {
    print "Usage: atacshift.pl input.sam genome.fa.sizes\n";
    exit;
}

# put the genome.fa.sizes in a hash
open my $in_sizes, "<:encoding(utf8)", $ARGV[1];
while (my $line = <$in_sizes>) {
  chomp $line;

  my @split_line = split '\t', $line;
  $sizes{$split_line[0]} = $split_line[1];
}

# declare variables
my @read;
my @cigar;
my $read_len;

# for each read...
open my $in_sam, "<:encoding(utf8)", $ARGV[0];
while (my $line = <$in_sam>) {
  chomp $line;

  # ignore header lines
  if ($line =~ (/^(\@)/)) {
    print "$line\n";
    next;
  }

  # split line into array
  @read = split '\t', $line;

  # only keep (necessary) first 11 columns
  @read = @read[0..10];
  
  if (not ($read[1] & 0x04)) {  # read must be mapped to shift
    # empty QUAL and SEQ
    $read[9] = "*";
    $read[10] = "*";

    @cigar = ( $read[5] =~ /\d+/g );
    my ($var1, $var2) = $read[5] =~ /(.*[A-RT-Z])(\d+S)?/; 
    $read_len = sum(split /[A-Z]/, $var1);
    if ($read[1] & 0x10) {  # reverse strand
      if (($read[1] & 0x01) && not ($read[1] & 0x08)) {  # read paired and mate not unmapped
          $read[7] = $read[7] + 4;  # mate shifts 4
          $read[8] = $read[8] + 9;  # template length decreases 4 + 5
      }
      $read_len += -5;
    }
    else {  # forward strand
      if (($read[1] & 0x01) && not ($read[1] & 0x08)) {  # read paired and mate not unmapped
          $read[8] = $read[8] - 9;  # template length decreases 4 + 5
      }
      $read[3] = $read[3] + 4;  # shift the position +4
      $read_len += -4;
    }
    $read[5] = "${read_len}M";

    # if the read falls outside of the chromosome after shifting; ignore it
    if (($read[3] < 0) || ($read[3] >= $sizes{$read[2]})) {
      next;
    }
  }
  print join "\t", @read, "\n";
}
