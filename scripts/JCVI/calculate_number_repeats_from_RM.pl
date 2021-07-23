#!/usr/local/bin/perl
use strict;
use Data::Dumper;

my $usage = "$0 -i repeat_masker_out_file -bin # [1-100 %, default 10] -id Repeat_id [defaul all] -gap allowed_gap [default 50]\n\n";
my %arg = @ARGV;

die $usage if $ARGV[0] eq '-h' or scalar @ARGV == 0;
my $gap = $arg{-gap} || 50;
my $rep_id = $arg{-id} || '.+';
my $binner = $arg{-bin} || 10;
my $file = $arg{-i} || die $usage;
my ($end_complete, $end_partial, $complete, $partial, $longest_hit, $shortest_hit, $per_coverage) = (0,0,0,0,0, 1000000,0);
my %bin_counter;
my $rep_len = 0;


open(FHI,"<$file") || die "ERROR, I cannot open $file: $!\n\n";

my @row;

while(<FHI>){
	chomp;
	s/^\s+//;
	s/\(|\)//g;
	my $row = get_row($_);	
	next unless $row->{rep_id} =~ m/$rep_id/;
	#if($row->{hit_len} > 10000){print "$_\n";exit};
	#next unless $row->{asmbl_id} == 2808 and $row->{rep_id} eq 'LINE1';
	push @row, $row;
}

my $old_row = shift @row;
my @processed_data;

## fuse overlapping hits
while (my $new_row = shift @row){
	my $comparisson = compare_hits($old_row, $new_row, $gap);
	if ($comparisson){
		## rows overlapp => $comparisson contains fused rows
		$old_row = $comparisson;
		push @processed_data, $old_row if scalar @row == 0; ## load last element of array
	}
	else {
		## rows do not overlap => load $old_row into @processed_data array
		push @processed_data, $old_row;
		push @processed_data, $new_row if scalar @row == 0; ## load last element of array
		$old_row = $new_row;
	}
}

## analyze processed hits
foreach my $row (@processed_data){
	
	#print Dumper($row) if $row->{hit_len} < 0;
	
	## Calculates longest & shortest hits
	($longest_hit, $per_coverage) = $row->{hit_len} >= $longest_hit ? ($row->{hit_len},$row->{hit_percent})  : ($longest_hit,$per_coverage);
	$shortest_hit = $row->{hit_len} <= $shortest_hit ? $row->{hit_len} : $shortest_hit;
	
	## Get longest repeat length
	$rep_len = $row->{rep_len} >= $rep_len? $row->{rep_len} : $rep_len;
	
	## Process bins
	my $bin = 0;
	do {
		$bin += $binner;
		print STDERR "ERROR, bin greater than 100%!! ($row->{hit_percent})\n\n" if $bin > 100 and $row->{hit_percent} > 100;
		
		$bin_counter{$bin}++ if $row->{hit_percent} <= $bin;
	} while($row->{hit_percent} >= $bin);
	
	## If repeat is @ the end of asmbl
	if($row->{asmbl_5} <= 5 and $row->{asmbl_left} <= 5){
		##################################
		# NOTE I don't count this because it is hard to tell if the asmbly is not a duplication
		##################################
		## if truncated at both end because end of asmbly then count as end_complete
		#$end_complete++ if $row->{hit_percent} < 90;
		#$complete++ if $row->{hit_percent} >= 90;
	}
	elsif($row->{asmbl_5} <= 5){
		## truncated at just 5' end of asmbly
		if( $row->{strand} eq '+' and $row->{rep_left} <= 5){
			$end_complete++ if $row->{hit_percent} < 90;
			$complete++ if $row->{hit_percent} >= 90;		
		}
		elsif( $row->{strand} eq 'C' and $row->{rep_5} <= 5){
			$end_complete++ if $row->{hit_percent} < 90;
			$complete++ if $row->{hit_percent} >= 90;					
		}
		else {
			$end_partial++;
		}
	}
	elsif($row->{asmbl_left} <= 5){
		## truncated at just 3' end of asmbly
		if( $row->{strand} eq '+' and $row->{rep_5} <= 5){
			$end_complete++ if $row->{hit_percent} < 90;
			$complete++ if $row->{hit_percent} >= 90;					
		}
		elsif( $row->{strand} eq 'C' and $row->{rep_left} <= 5){
			$end_complete++ if $row->{hit_percent} < 90;
			$complete++ if $row->{hit_percent} >= 90;					
		}
		else {
			$end_partial++;
		}	
	}
	## hit is internal
	elsif($row->{hit_percent} >= 90){
		$complete++;
	}
	else {
		$partial++;
	}
}
close FHI;
my $total_complete = $complete;
my $total_partial = $end_complete + $end_partial + $partial;
my $total = $total_complete + $total_partial;

print "$rep_id statistics:\n\n";
print "Number of complete truncated by end of asmbly = $end_complete\n";
print "Number of partial and truncated by end of asmbly = $end_partial\n";
print "Number of complete elements (>= 90% of the length) = $complete\n";
print "Number of partial elements = $partial\n\n";
print "Number of total complete elements = $total_complete\n";
print "Number of total partial elements (including complete truncated by end of assembly) = $total_partial\n";
print "Number of total elements = $total\n\n";
print "Longest hit = $longest_hit bp ($per_coverage % coverage)\n";
print "Shortest hit = $shortest_hit bp\n\n";
print "Bin:\n";
foreach my $bin (sort{ $a <=> $b} keys %bin_counter){
	my $bin_low = $bin - $binner;
	my $hit_len_low = int($rep_len*$bin_low/100);
	my $hit_len_top = int($rep_len*$bin/100);
	print "hits between $bin_low% ($hit_len_low bp) and $bin% ($hit_len_top bp) of coverage = ",$bin_counter{$bin},"\n";
}

######################################################################################################
sub compare_hits {
	my ($or, $nr, $gap) = @_;
	
	## possitive strand
	if( $or->{strand} eq $nr->{strand} 
			and $or->{strand} eq '+' 
			and $or->{r_mid_point} < $nr->{r_mid_point} 
			and $nr->{asmbl_5} - $or->{asmbl_3} <= $gap 
			and $or->{rep_id} eq $nr->{rep_id}
			and $or->{asmbl_id} eq $nr->{asmbl_id}
		){
		## Overlapped/contiguous hits
		my ($new_asmbl_end5, $ax, $ay, $new_asmbl_end3) = sort {$a <=> $b} ($or->{asmbl_5}, $or->{asmbl_3}, $nr->{asmbl_5}, $nr->{asmbl_3});
		my ($new_rep_end5, $rx, $ry, $new_rep_end3) = sort {$a <=> $b} ($or->{rep_5}, $or->{rep_3}, $nr->{rep_5}, $nr->{rep_3});
		my $asmbl_mid_point = int( ( $new_asmbl_end5 + $new_asmbl_end3 ) /2);
		my $rep_mid_point = int( ( $new_rep_end5 + $new_rep_end3 ) /2);
		
		my ($overlap,$x) = sort {$b<=>$a} ( ($or->{rep_3} - $nr->{rep_5}),($or->{asmbl_3} - $nr->{asmbl_5}) ); 
		$overlap = 0 if $overlap < 0;
		
		my $hit_len = $or->{hit_len} + $nr->{hit_len} - $overlap;
		my $overlap_percent = int( $overlap * 100 / $or->{rep_len});
		my $hit_percent =  $or->{hit_percent} + $nr->{hit_percent} - $overlap_percent;
		
		return {	asmbl_id => $or->{asmbl_id},
					rep_id => $or->{rep_id},
					asmbl_5 => $new_asmbl_end5,
					asmbl_3 => $new_asmbl_end3,
					asmbl_left => $nr->{asmbl_left},
					strand => $or->{strand},
					rep_5 => $new_rep_end5,
					rep_3 => $new_rep_end3,
					rep_left => $nr->{rep_left},
					rep_len => $or->{rep_len},
					hit_percent => $hit_percent,
					hit_len => $hit_len,
					a_mid_point => $asmbl_mid_point,
					r_mid_point => $rep_mid_point};
	}
	
	## negative strand
	elsif ( $or->{strand} eq $nr->{strand} 
			and $or->{strand} eq 'C' 
			and $or->{r_mid_point} > $nr->{r_mid_point} 
			and $nr->{asmbl_5} - $or->{asmbl_3} <= $gap 
			and $or->{rep_id} eq $nr->{rep_id}
			and $or->{asmbl_id} eq $nr->{asmbl_id}
		){
		## Overlapped/contiguous hits
		my ($new_asmbl_end5, $ax, $ay, $new_asmbl_end3) = sort {$a <=> $b} ($or->{asmbl_5}, $or->{asmbl_3}, $nr->{asmbl_5}, $nr->{asmbl_3});
		my ($new_rep_end5, $rx, $ry, $new_rep_end3) = sort {$a <=> $b} ($or->{rep_5}, $or->{rep_3}, $nr->{rep_5}, $nr->{rep_3});
		my $asmbl_mid_point = int( ( $new_asmbl_end5 + $new_asmbl_end3 ) /2);
		my $rep_mid_point = int( ( $new_rep_end5 + $new_rep_end3 ) /2);
		
		my ($overlap,$x) = sort {$b<=>$a} ( ($nr->{rep_3} - $or->{rep_5}),($or->{asmbl_3} - $nr->{asmbl_5}) ); 
		$overlap = 0 if $overlap < 0;
		my $hit_len = $or->{hit_len} + $nr->{hit_len} - $overlap;
		my $overlap_percent = int( $overlap * 100 / $or->{rep_len});
		my $hit_percent =  $or->{hit_percent} + $nr->{hit_percent} - $overlap_percent;
		
		
		return {	asmbl_id => $or->{asmbl_id},
					rep_id => $or->{rep_id},
					asmbl_5 => $new_asmbl_end5,
					asmbl_3 => $new_asmbl_end3,
					asmbl_left => $nr->{asmbl_left},
					strand => $or->{strand},
					rep_5 => $new_rep_end5,
					rep_3 => $new_rep_end3,
					rep_left => $nr->{rep_left},
					rep_len => $or->{rep_len},
					hit_percent => $hit_percent,
					hit_len => $hit_len,
					a_mid_point => $asmbl_mid_point,
					r_mid_point => $rep_mid_point};
	
	}
	else {
		return 0;
	}
		
	## negative strand
}

sub get_row {
	my @row = split(m/[\s\t]+/,$_[0]);
	my $rep_left = $row[8] eq '+' ? $row[12] : $row[10];
	
	if($row[8] eq '+'){
		my $hit_len = abs($row[10] - $row[11]) + 1;
		my $rep_len = $row[11] + $row[12];
		my $hit_percent =  int($hit_len * 100 / $rep_len);
		my $asmbl_mid_point = int( ($row[5] + $row[6] ) /2);
		my $rep_mid_point = int( ($row[10] + $row[11] ) /2);
		return {	asmbl_id => $row[4],
					rep_id => $row[9],
					asmbl_5 => $row[5],
					asmbl_3 => $row[6],
					asmbl_left => $row[7],
					strand => $row[8],
					rep_5 => $row[10],
					rep_3 => $row[11],
					rep_left => $rep_left,
					rep_len => $rep_len,
					hit_percent => $hit_percent,
					hit_len => $hit_len,
					a_mid_point => $asmbl_mid_point,
					r_mid_point => $rep_mid_point};
	} 
	elsif ($row[8] eq 'C') {
		my $hit_len = abs($row[12] - $row[11]) + 1;
		my $rep_len = $row[11] + $row[10];
		my $hit_percent =  substr($hit_len * 100 / $rep_len,0,4);
		my $asmbl_mid_point = int( ($row[5] + $row[6] ) /2);
		my $rep_mid_point = int( ($row[12] + $row[11] ) /2);
		return {	asmbl_id => $row[4],
					rep_id => $row[9],
					asmbl_5 => $row[5],
					asmbl_3 => $row[6],
					asmbl_left => $row[7],
					strand => $row[8],
					rep_5 => $row[12],
					rep_3 => $row[11],
					rep_left => $rep_left,
					rep_len => $rep_len,
					hit_percent => $hit_percent,
					hit_len => $hit_len,
					a_mid_point => $asmbl_mid_point,
					r_mid_point => $rep_mid_point};
	}
}
