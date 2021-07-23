#!/usr/local/bin/perl

use strict;


my $usage = "$0 < multiFastaFile > STDOUT";
my %arg = @ARGV;
die $usage if exists $arg{-h};
my $seq_len = 0;
my $seq;
my @array = qw/A C D E F G H I K L M N P Q R S T V W Y/;
my %amino;
my %percent;

while(<STDIN>){
	chomp;
	if (m/^>/){
		my $id = $_ ;
	} else {
		$seq = $_;
		$seq_len += length $seq;
		foreach my $aa (@array){
			while($seq =~ /$aa/ig){
				$amino{$aa}++;
			}
		}
	}
}

## calc %
foreach my $aa (@array){
	$percent{$aa} = $amino{$aa}/$seq_len*100;
	$percent{$aa} =~ m/(\d+\.\d{2})/;
	$percent{$aa} = $1;
}
foreach my$aa (@array){
	print "$aa: $amino{$aa} $percent{$aa}%\n";
}
