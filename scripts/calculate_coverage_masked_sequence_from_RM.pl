#!/usr/local/bin/perl
use strict;

my $usage = "$0 -r <RepeatMasker out file> -a <asmbl asmbl_len tab file>\n\n";
my %arg = @ARGV;
die $usage if ( scalar(@ARGV) == 0 || $ARGV[0] eq '-h' || !$arg{-r} || !$arg{-a});

open(RM,"<$arg{-r}") || die "ERROR, I cannot open $arg{-r}: $!\n\n";

my %asmbl_list;
my %repeat_list;
my %coverage;

while(<RM>){
	chomp;
	next if m/^[\s\t]*$/;
	s/^\s+//;
	my $line = &get_repeat_info($_);
	$asmbl_list{ $line->{asmbl} } = $line->{asmbl_len};
	$repeat_list{ $line->{repeat} } = 1;
	
	# initialize %coverage
	$coverage{ $line->{asmbl} }->{ $line->{repeat} } = 0 unless $coverage{ $line->{asmbl} }->{ $line->{repeat} };
	
	
	$coverage{ $line->{asmbl} }->{ $line->{repeat} } += $line->{coverage};
	$coverage{ $line->{asmbl} }->{total} += $line->{coverage};
	
}
close RM;

open(ASMBL,"<$arg{-a}") || die "ERROR, I cannot open $arg{-a}: $!\n\n";
while(<ASMBL>){
	chomp;
	next if m/^[\s\t]*$/;
	my ($asmbl, $asmbl_len) = split m/\t+/;
	warn "WARNING, wrong asmbly file format ($asmbl, $asmbl_len)\n", next unless ($asmbl && $asmbl_len);
	
	if (exists $asmbl_list{$asmbl} ){
		## print coverage for each repeat
		foreach my $repeat (keys %repeat_list){
			next unless $repeat;
			if (exists $coverage{$asmbl}->{$repeat}){
				print "$asmbl\t$asmbl_len\t$coverage{$asmbl}->{$repeat}\t$repeat\n";
			}
		}
		## print total coverage
		print "$asmbl\t$asmbl_len\t$coverage{$asmbl}->{total}\ttotal\n";
	}
	else {
		print "$asmbl\t$asmbl_len\t0\ttotal\n";
	}
	print "\n";
}
close ASMBL;

############################################################

sub get_repeat_info {
	my $line = shift;
	my @line = split(m/\s+/,$line);
	my $asmbl_left = $1 if $line[7] =~ m/(\d+)/;
	my $asmbl_len = $asmbl_left + $line[6];
	my $hit_len = $line[6] + 1 - $line[5];
	my $coverage = $1 if ($hit_len * 100 / $asmbl_len) =~ m/^(\d+)\./;
	return {	asmbl => $line[4],
				asmbl_len => $asmbl_len,
				repeat => $line[9],
				strand => $line[8],
				hit_len => $hit_len,
				end5 => $line[5],
				end3 => $line[6],
				coverage => $coverage,
				};
}
