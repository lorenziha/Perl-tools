#!/Users/lorenziha/bin/miniconda3/bin/perl
use strict;

my (%h);
while(<>){
	chomp;
	next if m/^[\s\t]*$/;
	$h{$_}++;
}
foreach my $key (keys %h){
	print "$key\t$h{$key}\n";
}
