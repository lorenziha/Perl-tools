#!/usr/bin/perl
use strict;

my $usage = "$0 -i <allele file> [-c <minimum read coverage of allele [2]> -aa <print only alt alleles T/F [F]> -q <mean read quality [10]> -rb <min fwd/rev read ratio [0.1]> -rt <max fwd/rev read ratio [0.9]>]\n\n";

my %arg = @ARGV;
die $usage unless $arg{-i};	

my $QUAL = $arg{-q} || 10;
my $RATIO_BOTTOM = $arg{-rb} || 0.1;
my $RATIO_TOP = $arg{-rt} || 90.;
my $AA_FLAG = $arg{-aa} eq 'T'? 1 : 0;
my $MIN_COV = $arg{-c} || 2;

print "#Chromosome\tPosition\tA\tT\tC\tG\n";

my ($chrom, $pos, $allele_info, %hist );
open(ALLELE, "<$arg{-i}") || die "ERROR, I cannot open $arg{-i}: $!\n\n";
while(<ALLELE>){
	chomp;
	next if m/^[\s\t]*$/;
	unless(m/^\d/){
		$chrom = $_;
		next;
	}
	my ($pos, $allele_info) = split /\t/;
	my ($a, $t, $c, $g, $fwd, $rev, $qual) = split(";", $allele_info);
	my $total = $a + $t + $c + $g;
	next if $total == 0; # skip positions with zero reads.
	
	# Filter out allele supported for less than $MIN_COV
	next if (($a>0 && $a < $MIN_COV) || ($t >0 && $t < $MIN_COV) || ($c >0 && $c < $MIN_COV) || ($g >0 && $g < $MIN_COV) );
	
	# calculate allele freqs
	my $fa = sprintf("%.0f", 100 * ($a/$total));
	my $ft = sprintf("%.0f", 100 * ($t/$total));
	my $fc = sprintf("%.0f", 100 * ($c/$total));
	my $fg = sprintf("%.0f", 100 * ($g/$total));
	my $max_freq = 	&max($fa, $ft, $fc, $fg);

	# calculate average read qualitymy 
	my $av_qual = $qual/$total;

        # calculate forward/reverse ratio
        my $fr = $fwd/$total;
        
       
	
	# Calculate histogram of frequences
	foreach my $freq ( ($fa, $ft, $fc, $fg) ){
		$hist{$freq}++ unless ($freq == 0 || $freq == 100 || $av_qual < $QUAL); # ignore frequences for monoallelic positions.
	}

	# print "($fr > $RATIO_BOTTOM && $fr < $RATIO_TOP && $av_qual >= $QUAL)\n";

	if($fr >= $RATIO_BOTTOM && $fr <= $RATIO_TOP && $av_qual >= $QUAL){
		if($AA_FLAG){
			print "$chrom\t$pos\t$fa\t$ft\t$fc\t$fg\n" if($max_freq < 100); 
		} 
		else {
			print "$chrom\t$pos\t$fa\t$ft\t$fc\t$fg\n";
		}
	}
}
close ALLELE;

# Print out histogram of frequences
print "\n# Histogram\n\n";
foreach my $freq (sort {$a <=> $b} keys %hist){
	print "$freq\t$hist{$freq}\n";
}

sub max {
	my @x = @_;
	my @sort = sort {$a <=> $b} @x;
	#print "@sort\n";
	return $sort[3];
}
