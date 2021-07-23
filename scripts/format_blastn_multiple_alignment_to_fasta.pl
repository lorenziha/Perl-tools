#!/usr/local/bin/perl
use strict;

my $usage = "$0 -i blastn_multiple_alignment_output_file -d database_fasta_file\n\nNOTE: IDs have to start with a number!!\n\n";
my %arg = @ARGV;
die $usage unless $arg{-i} and $arg{-d};

open (FHI,"<$arg{-d}") || die "ERROR, I cannot open $arg{-d} : $!\n\n";
my @ID;
while(<FHI>){
      chomp;
      if(m/^>(\S+)/){
           push @ID, $1;
      }
}
close FHI;

open (FHI,"<$arg{-i}") || die "ERROR, I cannot open $arg{-i} : $!\n\n";
my (%h, %z);
while(<FHI>){
	chomp;
	next if m/^[\s\t]*$/;
	next unless m/^\d+/; 
	my @x=split /\s+/; 
	my $L=60-length($x[2]);
	my $s= "-" x $L.$x[2]; 
	$h{$x[0]} .= uc($s); 
	push @{$z{$x[0]}},$x[1],$x[3];  
} 

my $id;
foreach my $k (keys %h){
        if($k eq '1_0'){
               $id = "query";
        }
        else{
               $id = @ID[$k];
               #my $id="seq_$k $z{$k}[0]-$z{$k}[-1]" ; 
        }
        print STDOUT ">$id\n$h{$k}\n" ;
}



