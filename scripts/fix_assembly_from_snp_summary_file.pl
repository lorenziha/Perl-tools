#!/usr/local/bin/perl

use strict;
use lib ($ENV{EUK_MODULES});
use Fasta_reader;


my $usage = "$0 -f <multifasta_file> -s SNP_summary_file\n\n";

die $usage unless $ARGV[0];
my %arg = @ARGV;
my @keywords;
my $fasta_file = $arg{-f} || die $usage;
my $snp = $arg{-s} || die $usage;

my %snp;
open(SNP_FILE,"<$snp") || die "$!:$snp\n\n";
while(<SNP_FILE>){
	chomp;
	my ($id, $pos, $ref, $alt, @rest) = split /\t/;
	push @{$snp{$id}}, "$pos:$ref:$alt";
	
}
close SNP_FILE;

my $fasta_reader = new Fasta_reader($fasta_file);
my %new_seq;
while( my $seqObj = $fasta_reader->next() ){
	my $header = $seqObj->get_header();
	my $seq = $seqObj->get_sequence();
	if ($snp{$header}){
	my $offset = 0;
		foreach my $snp_data (@{$snp{$header}}){
			my ($pos,$ref,$alt) = split(":", $snp_data);
			my $len = length($ref);
			substr($seq, $pos + $offset - 1, $len, $alt);
			$offset += ( length($alt)-$len);
		}	

	}
	$new_seq{$header} = $seq;	

}

foreach my $id (keys %new_seq){
	print ">new_$id\n$new_seq{$id}\n";
}

