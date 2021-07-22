#!/Users/lorenziha/bin/miniconda3/bin/perl
use strict;

my $usage = "$0 -i <depth_file> -b <bin size in bp [120]> -r <rename toxo ctg IDs T/F [F]>\n\n";
my %arg = @ARGV;

die $usage unless $arg{-i};
my $BIN = $arg{-b} || 120;
my $RENAME = $arg{-r} eq 'T'? 1:0;

open (DEPTH, "<$arg{-i}") || die "ERROR, I cannot open $arg{-i}: $!\n\n";
my ($chr, $c, $cov, $coord);
while(<DEPTH>){
	chomp;
	my @x = split /\t/;
	
	$x[0] =~ s/^TGME49_//; $x[0]=~s/KE138841/apicoplast/ if $RENAME;

	if($chr ne $x[0]){
		$c = 0; 
		$cov = 0; 
		$coord = 0
	};
	$chr = $x[0];
	$c++; 
	$coord += $x[1];
	$cov += $x[2];
	if($c == $BIN){
		my $mean_cov = int($cov/$BIN); 
		if($mean_cov == 0){
			$mean_cov=0.00001
		};
		my $mean_coord = int($coord/$BIN); 
		print "$x[0]\t$mean_coord\t$mean_cov\n"; 
		$c = 0; 
		$cov = 0; 
		$coord = 0 
	}
}
close DEPTH;
exit(0);
