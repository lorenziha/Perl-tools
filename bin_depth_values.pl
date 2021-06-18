#!/usr/bin/perl
use strict;
use Statistics::Descriptive; 
my $usage = "$0 -d depth_file [-b bin_size_in_bp [1000]]\n\nDepth file format:\n\n   Chr_ID<TAB>Coord<TAB>Depth\n\n";
my %arg = @ARGV;

die $usage unless $arg{-d};

my $BIN = $arg{-b} || 1000;

open (DEPTH, "<$arg{-d}")|| die "ERROR, I cannot open $arg{-d}: $!\n\n";
my ($c, $chr, @coord, @cov, $mean_cov, $mean_coord);
while(<DEPTH>){
	chomp;
	my @x = split /\t/; 
	if($chr ne $x[0]){
		$c = 0; 
		@cov = (); 
		@coord = ();
	};
	$chr = $x[0];
	$c++; 
	push @coord, $x[1]; 
	push @cov,   $x[2];
	if($c == $BIN){ 
		my $stat_coord = Statistics::Descriptive::Full->new();
		$stat_coord->add_data(@coord);

		my $stat_cov = Statistics::Descriptive::Full->new();
		$stat_cov->add_data(@cov);

		$mean_cov = $stat_cov->median(); 
		$mean_coord = $stat_coord->mean(); 
		print "$x[0]\t$mean_coord\t$mean_cov\n"; 
		$c = 0; 
		@cov = (); 
		@coord = (); 
	}
}

close DEPTH;
