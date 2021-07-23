#!/usr/local/bin/perl

use strict;
use lib ($ENV{EUK_MODULES});
use Fasta_reader;

my $usage = "$0 -i multiFastaFile > STDOUT";
my %arg = @ARGV;
die $usage if exists $arg{-h};

my $fasta_file = new Fasta_reader($arg{-i}) || die "I cannot find input file\n";
my $count = 0;



my @array = qw/A C D E F G H I K L M N P Q R S T V W Y/;
my $fasta_obj;

print "aminoacids:\t\t",join("\t",@array),"\n";

while(my $fasta_obj = $fasta_file->next() ){
	chomp;
	my %amino;
	my %percent;
	my $id = $fasta_obj->get_accession();
	my $seq = $fasta_obj->get_sequence();
	my $len = length($seq);
	next if $len == 0;
	$count++;
	foreach my $aa (@array){
		$amino{$aa} = countString($seq,$aa);
	}
	
	## calc %
	foreach my $aa (@array){
		$percent{$aa} = ($amino{$aa}/$len)*100;
		if( $percent{$aa} =~ m/(\d+\.\d{2})/){
			$percent{$aa} = $1;
		}
	}


	print "$count\t$id:\t\t";
	foreach my $aa (@array){
		print "$amino{$aa}\t";
	}
	print "\n$count\t$id\t$len\t";
	foreach my $aa (@array){
		print "$percent{$aa}\t";
	}
	print "\n";

}

sub countString{
	my ($seq, $aa) = @_;
	my $counter = 0;
	while($seq =~ /$aa/ig){
		$counter++;
	}
	return $counter;
}
