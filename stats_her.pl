#!/usr/bin/perl
use strict;

use Statistics::Descriptive;

my $usage = "$0  <STDIN with list of numbers>\n\n";
my @data;
while(<>){
	chomp;
	next unless m/^[\s\t]*\d+\.?\d*/;
	push @data, $_ ;
}

my $stat = Statistics::Descriptive::Full->new();
$stat->add_data(@data);
my $mean = $stat->mean();
my $var  = $stat->variance();
my $count = $stat->count();
 
my $min = $stat->min();
my $max = $stat->max();
my $sum = $stat->sum();
my $sd  = $stat->standard_deviation();
my $median = $stat->median();
my $mode = $stat->mode();

print "==================================\nResults:\nCount = $count\nMean = $mean\nSD = $sd\nVar = $var\nSum = $sum\nMin = $min\nMax = $max\nMedian = $median\nMode = $mode\n\n";


