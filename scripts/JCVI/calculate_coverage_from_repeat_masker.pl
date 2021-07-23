#!/usr/local/bin/perl
use strict;
use Data::Dumper;

my $usage = "$0 -i repeat_masker_output -id repeat_id -type type_id [all = *] \n\n";

die $usage if $ARGV[0] eq '-h' or scalar @ARGV == 0;
my %arg = @ARGV;
my $file = $arg{-i} || die $usage;
my $rep_id = $arg{-id} || '.+';
my $type = 0;
if(exists $arg{-type}){
	$rep_id = $arg{-type};
	$type = 1;
}
my %asmbls;

open(FHI,"<$file") || die "ERROR, I cannot open $file: $!\n\n";
while(<FHI>){
	chomp;
	s/^\s+//;
	my @x = split m/[\s\t]+/;
	next unless $x[9+$type] =~ m/$rep_id/;
	push @{ $asmbls{ $x[4] } },{ end5 => $x[5], end3 => $x[6] };
}
close FHI;

my $coverage = 0;
foreach my $asmbl ( keys %asmbls){
	my @sorted = sort { $a->{end5} <=> $b->{end5} } (@{ $asmbls{$asmbl} });
	#print Dumper(@sorted);	
	if(scalar @sorted > 1){
		my $old_coords = shift @sorted;
		## Add first pair of coords per asmbly
		$coverage += ($old_coords->{end3} - $old_coords->{end5} + 1 );
		foreach my $coords (@sorted){
			my $new_coverage = 0;
			if($coords->{end5} <= $old_coords->{end3}){
				## overlap
				$new_coverage = abs($coords->{end3} - $old_coords->{end3});
				
			}
			else {
				## No overlap
				$new_coverage = $coords->{end3} - $coords->{end5} + 1;
			}
			$coverage += $new_coverage;
		}
	}
	elsif(scalar @sorted == 1){
		my $new_coverage = $sorted[0]->{end3} - $sorted[0]->{end5} + 1;
		$coverage += $new_coverage;
	}
	else {
		print STDERR "ERROR, no coordinates for $asmbl!\n\n";
	}
}

print "Total repeat coverage for $rep_id = $coverage bp\n\n";
