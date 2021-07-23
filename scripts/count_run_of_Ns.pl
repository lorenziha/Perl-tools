#!/usr/local/bin/perl

use strict;
use warnings;
use lib ($ENV{EUK_MODULES});
use Fasta_reader;

my $usage = "usage: $0 -f multifasta_file\n\n";

die $usage if $ARGV[0] eq '-h';

my %arg = @ARGV;
my $file = $arg{-f} || die $usage;

my $fasta_reader = new Fasta_reader($file);
my $seq_counter  = 0;
my $count_by_seq = 0;
my $count_total  = 0;
my $len_per_seq  = 0;
my $len_total    = 0;

while (my $seqObj = $fasta_reader->next()) {
	$seq_counter++;
	
    my $header  = $seqObj->get_header();
    my $seq     = $seqObj->get_sequence();

	while( $seq =~ m/(N+)/ig){ 
		$count_by_seq++; 
		$len_per_seq += length $1;
	}
	
	$count_total += $count_by_seq;
	$len_total   += $len_per_seq;
	print "$header : $count_by_seq ($len_per_seq bp)\n";
	$count_by_seq = 0;
	$len_per_seq  = 0;
}

print "\nTotal number of 'N' runs in $seq_counter seqs: $count_total ($len_total bp)\n";
exit(0);
