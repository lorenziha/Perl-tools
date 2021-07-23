#!/usr/local/bin/perl
use strict;

my $usage = "$0 -i <vcf file> -f min_allele_freq [0.9]>  -m min_alt_allele_counts [5] -c min_allele_coverage [10] -P <ploidy [1]> \n\n";
my %arg = @ARGV;
die $usage unless ( $arg{-i} );
my %names;
my $ALLELE_FREQUENCY = $arg{-f} || 0.9; ## minimum allele frequency to call variants in each sample
my $MIN_COV = $arg{-c} || 10;
my $ALT = $arg{-m} || 5;
my $PLOIDY = $arg{-P} || 1; ## default is haploid
die "ERROR, ploidy (-P) has to be either 1 or 2 (currently $PLOIDY)\n\n" unless $PLOIDY == 1 or $PLOIDY == 2;


open (FHI, "<$arg{-i}") || die "ERROR, I cannot open $arg{-i}: $!\n\n";
#print "CHROMOSOME\tPOSITION\tREFERENCE\tVARIANT\tCOVERAGE\tALLELE_FREQ\tQUAL_FILTER\tMUTATION_TYPE\tCDS_MUTATION\tPROT_MUTATION\tGENE_ID\tPRODUCT_NAME\n";
my (@freq, @depth);
while(<FHI>){
	if ( m/^#/){
		print "$_";
		next;
	}
	my $line = $_;
	chomp;
        my @x = split /\t/;
        my $depth_freq;
	my $COUNTER;
	my $ALT_COUNT;
	my ($chrom, $pos, $ref, $alt) = ($x[0], $x[1],$x[3], $x[4]);
	my @alleles = $ref;
	push @alleles, split("," , $alt);
	## Number of columns?
	my $l = scalar(@x)-1;
	foreach my $sample (@x[9..$l]){
	 
		my ($freq, $depth, $genotype) = &get_freq($sample);
		my $allele;
		if($PLOIDY == 1){
			#$allele = $freq >= $ALLELE_FREQUENCY ? $alt : $ref;
			$allele = @alleles[$genotype]; #  == 1? $alt : $ref;
		}
		if($PLOIDY == 2){
			my @predicted_alleles = split("/", $genotype);
			$allele = $alleles[$predicted_alleles[0]]."_".$alleles[$predicted_alleles[1]] ; #$freq <= $MIN ? "0_0" : $freq > $MAX ? "1_1" : "0_1";
		}
		#$depth_freq .= $depth."\t".$freq."\t$allele\t";
		next if $depth < $MIN_COV;
		if ($freq < $ALLELE_FREQUENCY ){
			$COUNTER;
		}
		elsif ($genotype > 0  ) {
			$ALT_COUNT++;
		}
	}
	print "$line" if ($COUNTER <= 4 && $ALT_COUNT >= $ALT) ;
}
close FHI;
exit(0);

##############################
sub get_freq {
	# GT:AD:DP:GQ:PL
	##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
	##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
	##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
	##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
	##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">

	my $data = shift;
	my @x = split(':',$data);
	
	my $AD = $x[1];
	my @AD =  split(',' , $AD); # my ($AD_R, $AD_A) = split(',' , $AD);
	my $DP = $x[2];
	my $GT = $x[0];
	#my ($r,$a)= ($x[2], $x[3]); # ($1, $2) if $x[1]=~m/^(\d+),(\d+)/;
	my $total = $DP ;
	#print "$data\n" if $total == 0;;
	my $freq;
	if ($total && $GT =~ m/^\d+$/){  
		$freq = $AD[$GT] / $DP ;
	} else {
		$freq = 0;
	}
	
	#print "a = $a , r = $r , $data\n";
	$freq = $1 if $freq =~ m/^(\d+\.?\d?\d?)/;
	return ($freq, $DP, $GT);
}

sub get_mutation {
	my $data = shift;
	my @stats = split(';', $data);
	my @annot = split(m/\|/, $data);
	#print "$stats[-1]\n";
	my ($pub_locus, $nuc_mutation, $mutation, $type) = ($annot[3], $annot[9],  $annot[10], $annot[1]);
	$pub_locus =~ s/\s+//g;
	my $my_mutation = $mutation || $nuc_mutation;
	return ($pub_locus, $nuc_mutation, $my_mutation, $type) ;
}

# INDEL;IDV=24;IMF=0.648649;DP=121;VDB=4.94589e-05;SGB=-2.73308;MQSB=0.120735;MQ0F=0;DPR=0,106;AF1=1;AC1=8;DP4=0,0,13,93;MQ=17;FQ=-90.9422;ANN=TGGG|frameshift_variant|HIGH|TGGT1_220650|TGGT1_220650|transcript|rna_TGGT1_220650-1|Coding|5/5|c.1033dupG|p.Glu345fs|1034/1110|1034/1110|345/369||WARNING_TRANSCRIPT_NO_START_CODON;LOF=(TGGT1_220650|TGGT1_220650|1|1.00)

