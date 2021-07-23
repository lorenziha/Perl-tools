#!/usr/local/bin/perl

use strict;
use lib ($ENV{EUK_MODULES});
use Fasta_reader;

## It searches a run of the specified aminoacid in a multi fasta file of protein sequences

my $usage = "$0 -i multiFastaFile -a Aminoacid_letter -w [windowsize|5] -n [number_of_ocurrencies|4] -x1 [relax_factor:0.1-1|1] > STDOUT";
my %arg = @ARGV;
die $usage if exists $arg{-h};

my $fasta_file = new Fasta_reader($arg{-i}) || die "I cannot find input file\n";
my @array = qw/A C D E F G H I K L M N P Q R S T V W Y/;
my $fasta_obj;
my $aa = $arg{-a} || die $usage;
my $window = $arg{-w} || 5;
my $ocurrencies = $arg{-n} || 4;
my $ratio = $ocurrencies/$window;
my $x1 = $arg{-x} || 1; ## relaxes the extension of the original seed changing the ratio by a factor = x1

die "wrong aa code!\n" unless grep(/$aa/,@array);
print "Seq_id\tStarting_position\tAA\tNum_od_Hits\tTotal_length\tSeq_length\tSequence\tRatio\tThreshold_ratio\n";
while(my $fasta_obj = $fasta_file->next() ){
	chomp;
	my $pos = 0;
	my $id = $fasta_obj->get_accession();
	my $seq = $fasta_obj->get_sequence();
	my $len = length($seq);
	next if $len < $window;
	
	while($pos < $len-$window){
		my $seed = substr($seq,$pos,$window);
		my $number_of_hits = countString($seed,$aa);
		my $seed_ratio = $number_of_hits/length($seed);
		if ($seed_ratio >= $ratio){
			my $counter = 0;
			while($seed_ratio >= ($ratio * $x1) ){
				$counter++;
				$seed = substr($seq,$pos,$window+$counter);
				$number_of_hits = countString($seed,$aa);
				$seed_ratio = $number_of_hits/length($seed);
				last if $pos+$window+$counter >= $len;
			}
			$seed = substr($seq,$pos,$window+$counter-1);
			$number_of_hits = countString($seed,$aa);
			$seed_ratio = $number_of_hits/length($seed);
			
			my $total_window = $window+$counter-1;
			print "$id\t$pos\t$aa\t$number_of_hits\t$total_window\t$len\t$seed\t$seed_ratio\t$ratio\n";
			$pos += $window+$counter-1;
			
		}
		$pos++;
	
	}
}

sub countString{
	my ($seq, $aa) = @_;
	my $counter = 0;
	while($seq =~ /$aa/ig){
		$counter++;
	}
	return $counter;
}
