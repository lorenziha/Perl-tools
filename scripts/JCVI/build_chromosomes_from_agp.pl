#!/usr/local/bin/perl
use strict;

my $usage = "$0 -i agp file -f contigs_file\n\n";
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

## read instruction file
my ($chrid, %chr, $scfseq);
open (FHI,"<$arg{-i}") || die "ERROR, I cannot open $arg{-i}:$!\n\n";
while(<FHI>){
	chomp;
	next if m/^[\s\t]*$/;
	if(m/^(chr\S+)/){
		$chrid = $1;
	}
	else {
		my @x = split /\t/;
		my $strand = $x[8];
		my $chrid = $x[0];
		my $typeid = $x[4];
		my $contigid = $x[5];
		if($typeid eq "N"){
			
			$scfseq = 'N' x $contigid;
		}
		else {
			$scfseq = $seq{$contigid};
		}
		if($strand eq '-'){
			$scfseq = reverse ($scfseq);
			$scfseq =~ tr/ACGTNacgtnyrkm/TGCANtgcanrymk/;
		}
		push @{$chr{$chrid}}, $scfseq;
	}
}

foreach my $k (keys %chr){
	print ">$k\n".join("", @{$chr{$k}})."\n";
}

sub get_seq {
	my ($id, $start, $end) = @_;
	$start = $start eq 'START'? 0 : $start-1;
	$end = $end eq 'END' ? length($seq{$id}) : $end;
	my $len = $end-$start;
	unless ($seq{$id}){die "ERROR, I cannot find sequence for ID \"$id\"\n\n"}
	return substr($seq{$id},$start,$len);	

}
close FHI;
