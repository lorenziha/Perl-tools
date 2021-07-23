#!/usr/local/bin/perl
use strict;

my $usage = "$0 -i wu-blast.btab [-M <Max distance between introns; default = 20000>]\n\n";
my %arg = @ARGV;
my $MAX_DIST_BETWEEN_HSPS = $arg{-M} || 20000;
die $usage unless exists($arg{-i}) ;

open(BTAB,"$arg{-i}") || die "ERROR, I cannot open $arg{-i}: $!\n\n";
my %query;
while(<BTAB>){
	chomp;
	my @x = split /\t/;
	push @{$query{$x[0].'-'.$x[5]}},[$x[2],$x[5], $x[6],$x[7],$x[8],$x[9],$x[10],$x[11],$x[19]];
}
close BTAB;

foreach my $query_asmbly (keys %query){
	my ($query, $asmbly) = split ("-",$query_asmbly);
	#print "$query, $asmbly, length query = $query{$query_asmbly}->[0][0]\n";
	my ($coverage);
	
	my $seq = 'o' x $query{$query_asmbly}->[0][0]; ## cretes artifactual sequence 
	#print "$seq\n";
	my @sorted_hsps = sort {$a->[4] <=> $b->[4]} @{$query{$query_asmbly}};
	my ($counter,$acc_old);
	my $identical;
	foreach my $hsp (@sorted_hsps){
		#print "$hsp->[4]\n";
		unless ($counter){$acc_old = $hsp->[4]};
		my $inter_hsp_dist = abs($acc_old - $hsp->[4]);
		if($inter_hsp_dist > $MAX_DIST_BETWEEN_HSPS){warn "distance between hits is too long\n\n"}	
		my $index = $hsp->[2]-1;
		my $len   = $hsp->[3] - $index;
		my $hit   = substr($seq,$index , $len);
		#print "before=$hit\n";
		$hit =~ s/H//g;
		#print "after=$hit\n";
		my $hit_len = length($hit); 
		substr($seq,$index , $len) = 'H' x $len;
		$identical += int($hsp->[7]*($hit_len)/100);
		#print "$hsp->[7] $identical\n";
	}
	#print "$seq\n";
	
	my $seqlen = length($seq);
	my $perc_id = ($identical / $seqlen)*100;
	#print "$perc_id\n";
	$perc_id =~ s/\.?\d+$//;
	$perc_id = $perc_id? $perc_id : 0;
	$seq =~ s/o+//g;
	my $shadowlen = length($seq);
	#print "$shadowlen / $seqlen\n";
	my $perc_cov = ($shadowlen/$seqlen)*100;
	$perc_cov =~ s/\.?\d+$//;
	$perc_cov = $perc_cov? $perc_cov:0;
	print "# query subject average_perc_identity percent_coverage\n";
	print "$query\t$asmbly\t$perc_id"."%\t$perc_cov"."%\n";
}



