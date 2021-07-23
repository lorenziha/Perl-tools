#!/usr/local/bin/perl


use lib ($ENV{EUK_MODULES}, $ENV{EGC_SCRIPTS});
use Egc_library;
use strict;
use DBI;
use Gene_obj;
use DBmodel_to_geneobj;
use Data::Dumper;

$|++;


my %asmbl;
my $family_id = '';
my $usage = "\n$0 -p <pub_locus file> -d <database> -c <chromosome maps> -s <scale_value> [0.01 - 1] -n <normailization_value> -l ruler >\n\n";

## Load Parameters

my %arg = @ARGV;
my $pub_locus_file = $arg{-p} || die $usage;
my $db = $arg{-d} || die $usage;
my $chromosomes = $arg{-c} || die $usage;
my $scale = $arg{-s} || 1;
if ($scale > 1 or $scale < 0.01){die $usage;}
my $normalize = $arg{-n} if exists ( $arg{-n} );

## Group pub_locus values by family_id from input file in HoA %asmbl 

open(PUB_LOCUS,"$pub_locus_file") || die "$usage $!\n";

my $family_count = 0;
my %family_number;

while(<PUB_LOCUS>){
	chomp;
	next unless m/[\d+\S+]/;
	if (m/^\s*Pv\d{6}/){
		s/\s+//g;
		die "Family ID is absent from input fiel\n" unless $family_id;
		push @{ $asmbl{$family_id} }, $_;
	}
	elsif (m/^\S+/){
		s/\s+//g;
		$family_id = $_;
		$family_count++;
		$family_number{$family_id} = $family_count;
	}
	else {
		die "Wrong input file format!\n";
	}
}
close PUB_LOCUS;

## Format chromosome asignment file and store in HoAoA %chromo

my %chromo;

open(CHROMO,"$chromosomes") || die "$usage $!\n";

while(<CHROMO>){
	chomp;
	if (m/^\S+\s(\d{1,3}):\s+(\S+)/){
		my $chromo_id = $1;
		my @path = split (m/\.{3,}/,$2);
		my $strand;
		my $asmbl_id;
		foreach (@path){
			if(m/B(\d+)E/){ ## positive strand
				$strand = '+';
				$asmbl_id = $1; 			
			}
			elsif (m/E(\d+)B/){ ## negative strand
				$strand = '-';
				$asmbl_id = $1;
			}
			elsif (m/(\d+)kb/i){
				$strand = $1;
				$asmbl_id = 'spacer';
			
			}
			else {next;}
			push @{ $chromo{$chromo_id} }, [$asmbl_id,$strand];
		}
		
	}
}
close CHROMO;

my %asmbl_len;
my %asmbl_genes;
my $flag;

foreach my $family (keys %asmbl){
	foreach my $pub_locus ( @{ $asmbl{$family} } ){
		
		## SQL QUERY: retrieve gene information
		
		my $query1 = "    SELECT distinct af.asmbl_id, af.end5
                	          FROM evidence e, clone_info ci, ident i, feat_link fl, 
                        	       asm_feature af, phys_ev pe
	                          WHERE af.feat_name = e.feat_name
				  AND fl.parent_feat = i.feat_name
				  AND fl.child_feat = af.feat_name
        	                  AND ci.asmbl_id = af.asmbl_id
                	          AND ci.asmbl_id = af.asmbl_id
                        	  AND pe.feat_name = af.feat_name
	                          AND pe.ev_type = 'working'
        	                  AND ci.is_public = 1
				  AND i.pub_locus = '$pub_locus'";
		
		my $query2 = "SELECT asmbl_id, length FROM clone_info WHERE is_public = 1";
			  
		my $dbh = DBI -> connect( "dbi:Sybase:server=SYBTIGR; packetSize=8092", 'access', 'access' );
		$dbh -> do( "use $db" );
		my $sth = $dbh -> prepare( $query1 );
	
        	## Do SQL QUERY
		
		my $asmbl_id;
		my $end5;
		my $asmbl_len;
		
		$sth -> execute();
		$sth -> bind_columns( \$asmbl_id, \$end5 );

		## Print query output

		while( $sth -> fetch() ){
	    		#$asmbl_len{$asmbl_id} = $asmbl_len;
			push @{ $asmbl_genes{$asmbl_id} }, [$pub_locus, $family, $end5];
		}
	
		$sth = $dbh -> prepare( $query2 );
	
        	## Get asmbl lengths

		if($flag < 1){		
			$sth -> execute();
			$sth -> bind_columns( \$asmbl_id, \$asmbl_len );

			## Print query output

			while( $sth -> fetch() ){
	    		$asmbl_len{$asmbl_id} = $asmbl_len;
			}		
		}
		$flag = 1;
		
		# Close DB
	
		$sth -> finish;
		$dbh -> disconnect;
	}

}

## Calculates total chromosome length

my %chromo_len;

foreach my $chromo_id (sort {$a <=> $b} keys %chromo){
	
	foreach my $array ( @{ $chromo{$chromo_id} } ){
		my $asmbl_id = $array->[0];
		my $orientation = $array->[1];
		if ($asmbl_id eq 'spacer'){
			$chromo_len{$chromo_id} += $orientation;
			next;
		}
		$chromo_len{$chromo_id} += $asmbl_len{$asmbl_id};
		#print "->$asmbl_id: $asmbl_len{$asmbl_id} ($chromo_len{$chromo_id})\n";
	}
}

## Calculates gene coordinates

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


#############################################################################
## begin program

