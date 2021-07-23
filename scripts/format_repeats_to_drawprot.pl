#!/usr/local/bin/perl


use lib ($ENV{EUK_MODULES}, $ENV{EGC_SCRIPTS});
use Egc_library;
use strict;
use DBI;
use Gene_obj;
use DBmodel_to_geneobj;
use Data::Dumper;

$|++;


my $usage = "\n$0 -r <RepeatMasker.out file> -i <repeat_id> -c <chromosome maps> -s <scale_value> [0.01 - 1] -n <normailization_value> -l ruler >\n\n";

## Load Parameters

my %arg = @ARGV;
die $usage if ($ARGV[0] eq '-h' || scalar(@ARGV) == 0);

my $RepeatMasker_file = $arg{-r} || die $usage;
my $chromosomes = $arg{-c};# || die $usage;
my $scale = $arg{-s} || 1;
#die $usage if ($scale > 1 or $scale < 0.01);
my $normalize = $arg{-n} if exists ( $arg{-n} );

## Group pub_locus values by family_id from input file in HoA %asmbl 

open(RM,"<$RepeatMasker_file") || die "I cannot open $RepeatMasker_file: $!\n";

my $repeat_count = 0;
my %repeat_number;
my %asmbl; #stores RM entries grouped by asmbl_id
my %repeat_id;

while(<RM>){
	chomp;
	next if m/^[\s\t]*$/;
	s/^\s+//;
	my $repeat_hit = &get_repeat_hit_info($_);
	$repeat_id{ $repeat_hit->{repeat} } = 0;
	push @{ $asmbl{$repeat_hit->{asmbl} } }, $repeat_hit;
}

# compute repeat_IDs
my $counter;
foreach my $repeat (keys %repeat_id){
	$counter++;
	$repeat_id{$repeat} = $counter;
}

foreach my $asmbl_id (sort{$a <=> $b} keys %asmbl ){
	&print_asmbly( @{ $asmbl{$asmbl_id} } );
} 
	

if(0){
	## Calculates gene coordinates
	my (%chromo,%asmbl_len, %asmbl_genes, %family_number,%chromo_len);
	foreach my $chromo_id (sort {$a <=> $b} keys %chromo){
		my $spacer = 0;
		my $len = $chromo_len{$chromo_id};
		$scale = $normalize/$chromo_len{$chromo_id} if exists($arg{-n});
		my $ruler = int(100000 * $scale); ## spacings in the ruler
		$chromo_len{$chromo_id} = int($scale * $chromo_len{$chromo_id} );

		print ">Chromosome:$chromo_id,$chromo_len{$chromo_id}\n";

		if( exists($arg{-l}) ){
			for (my $i=0; $i <= $chromo_len{$chromo_id}; $i += $ruler ){
				my $y = $i + 2;
				print "$i,$y,100kb_scale,20,r\n";
			}
		}
		foreach my $array ( @{ $chromo{$chromo_id} } ){
			my $asmbl_id = $array->[0];
			my $orientation = $array->[1];
			if ($asmbl_id eq 'spacer'){
				$spacer += $orientation;
				next;
			}
			my $asmbl_length = $asmbl_len{$asmbl_id};
			foreach my $gene ( @{ $asmbl_genes{$asmbl_id} } ){
				my $gene_id = $gene->[0];
				my $gene_family = $gene->[1];
				my $gene_pos = 0;
				if ($orientation eq '+'){
					$gene_pos = $spacer + $gene->[2];
				} 
				else {
					$gene_pos = $spacer + $asmbl_length - $gene->[2];
				}
				my $gene_pos2;
				$gene_pos = int( $scale * $gene_pos );
				if ( exists( $arg{-n} ) ){
					$gene_pos2 = $gene_pos + 5;
				} else {
					$gene_pos2 = int( $scale * 500) + $gene_pos;
				}
				#print "$chromo_id ($chromo_len{$chromo_id} bp), $asmbl_id, $gene_id, $gene_family, $gene_pos\n" 
				print "$gene_pos,$gene_pos2,$family_number{$gene_family},a,r,$gene_family\n"; 
			}
			$spacer += $asmbl_length;
		}
		print "\n";
	}
}

#############################################################################
## subroutines

sub scale_coords {
	my $hit = shift;
	$scale = $normalize/$hit->{a_len} if exists($arg{-n});
	my $ruler = int(100000 * $scale); ## spacings in the ruler
	$hit->{a_len} = int( $scale * $hit->{a_len} );
	$hit->{a_end_5} = int( $scale * $hit->{a_end_5} );
	$hit->{a_end_3} = int( $scale * $hit->{a_end_3} );
	$hit->{a_end_3}++ if $hit->{a_end_3} == $hit->{a_end_5};
	return $hit;
}

sub print_asmbly {
	my @hits = @_;
	
	#print asmbly header
	my $scaled_asmbly = &scale_coords($hits[0]);
	print ">$scaled_asmbly->{asmbl},$scaled_asmbly->{a_len}\n";
	
	
	
	#print entry line
	foreach my $hit (@hits){
		# scale coords
		my $scaled_hit = &scale_coords($hit);
		
		#print results
		print 	"$scaled_hit->{a_end_5},$scaled_hit->{a_end_3},".$repeat_id{ $scaled_hit->{repeat} }
				.",a,r,$scaled_hit->{repeat}\n";
	}
	
	## print a blanck line
	print "\n";
}

sub get_repeat_hit_info {
	my $hit = shift;
	my @line = split(m/\t+/, $hit);
	my $asmbl_leftover = $1 if $line[7] =~ m/(\d+)/;
	my $len = $asmbl_leftover + $line[6];
	my $strand = $line[8] eq '+' ? '+' : '-';
	return {
		asmbl => $line[4],
		a_end_5 =>$line[5],
		a_end_3 => $line[6],
		repeat => $line[9],
		a_len => $len,
		a_strand => $strand,
	};
}
