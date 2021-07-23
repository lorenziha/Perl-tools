#!/usr/local/bin/perl
use strict;

my $usage = "$0 -p <chromosome prefix>-i <output of map_asmbls_to_chromosomes_from_coords.pl> -f multifasta_file\n\n";
my %arg = @ARGV;
die $usage unless $arg{-i} and $arg{-f};

## Load scaffold sequences
my (%seq,$id,$seq);
open (FASTA,"<$arg{-f}") || die "ERROR, I cannot open $arg{-f}:$!\n\n";
while(<FASTA>){
	chomp;
	if(m/^>(\S+)/){
		$id = $1;
	}
	else {
		s/[\s\t]+//g;
		$seq{$id} .= $_;
	}

}
close FASTA; 

## read instruction file (chrID <TAB> queryID <TAB> strand(1/-1) )
## Instruction file is the output of map_asmbls_to_chromosomes_from_coords.pl
my ( %chr, %flag);
open (FHI,"<$arg{-i}") || die "ERROR, I cannot open $arg{-i}:$!\n\n";
while(<FHI>){
	chomp;
	next if m/^[\s\t]*$/;
	my ($ref, $query, $strand) = split m/\t/;
	## Replace prefix if defined
	if($arg{-p}){
		$ref =~ s/TGME49/$arg{-p}/;
	}
	if (!exists($flag{$ref})){
		## first scfld of chromosome
		$flag{$ref}++;
		print ">$ref\n";
		&print_seq($query,$strand);
	}
	else {
		print "N" x 100;
		&print_seq($query,$strand);
	}
}
exit(0);

##################################################
sub print_seq {
	my ($qid, $strand) = @_;
	if ($strand < 0){
		$seq{$qid} = reverse ($seq{$qid});
		$seq{$qid} =~ tr/ACGTNacgtnyrkm/TGCANtgcanrymk/;
	}
	print "$seq{$qid}\n";
}
