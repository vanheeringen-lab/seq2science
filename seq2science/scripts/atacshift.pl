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
my $rlen;

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

    # cut of trailing unmapped nucs
    @cigar = split /(?<=[A-Z])/, $read[5];
    @cigar = grep {!/\d(S|D)/}  @cigar; 
    #if ($cigar[0] =~ /\dS/) {
    #  shift @cigar;
    #}
    #if ($cigar[-1] =~ /\dS/) {
    #  pop @cigar;
    #}
    $rlen = sum(split /[A-Z]/, join '', @cigar);

    if ($read[1] & 0x10) {  # reverse strand
      if (($read[1] & 0x01) && not ($read[1] & 0x08) && not ($read[8] == 0)) {  # read paired and mate not unmapped
          $read[7] = $read[7] + 4;  # mate shifts 4
          $read[8] = $read[8] + 9;  # template length increases 4 + 5
      }
      $rlen -= 5;
      $read[9] = substr($read[9], 0, $rlen)
    }
    else {  # forward strand
      if (($read[1] & 0x01) && not ($read[1] & 0x08) && not ($read[8] == 0)) {  # read paired and mate not unmapped
          $read[8] = $read[8] - 9;  # template length decreases 4 + 5
      }
      $read[3] = $read[3] + 4;  # shift the position +4
      $rlen -= 4;
      $read[9] = substr($read[9], 4)
    }

    # make sure the starting position doesn't map off the chr
    if ($read[3] < 1) {
      $rlen -= 1 + $read[3];
      $read[3] = 1;
    }
    # make sure the cigar doesn't map off the chr
    if ($read[3] + $rlen >= $sizes{$read[2]}) {
      $rlen = $sizes{$read[2]} - $read[3] - 1;
    }

    $read[9] = substr($read[9], 0, $rlen);
    # Set the CIGAR
    $read[5] = "${rlen}M";

    # empty QUAL and SEQ
    # $read[9] = "*";
    $read[10] = "*";
  }
  print join "\t", @read, "\n";
}
