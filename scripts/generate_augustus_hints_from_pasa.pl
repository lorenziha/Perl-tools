#!/usr/local/bin/perl
use strict;

## INPUT

#373     alignAssembly-ME49_genes_tgrh88_pasa    cDNA_match      156928  157086  .       -       .       ID=chain_31;Target=asmbl_31 1558 1716 +
#373     alignAssembly-ME49_genes_tgrh88_pasa    cDNA_match      157320  157534  .       -       .       ID=chain_31;Target=asmbl_31 1343 1557 +
#373     alignAssembly-ME49_genes_tgrh88_pasa    cDNA_match      157641  158124  .       -       .       ID=chain_31;Target=asmbl_31 859 1342 +
#373     alignAssembly-ME49_genes_tgrh88_pasa    cDNA_match      158457  158666  .       -       .       ID=chain_31;Target=asmbl_31 649 858 +
#373     alignAssembly-ME49_genes_tgrh88_pasa    cDNA_match      159136  159458  .       -       .       ID=chain_31;Target=asmbl_31 326 648 +
#373     alignAssembly-ME49_genes_tgrh88_pasa    cDNA_match      160131  160256  .       -       .       ID=chain_31;Target=asmbl_31 200 325 +
#373     alignAssembly-ME49_genes_tgrh88_pasa    cDNA_match      160591  160695  .       -       .       ID=chain_31;Target=asmbl_31 95 199 +
#373     alignAssembly-ME49_genes_tgrh88_pasa    cDNA_match      161356  161414  .       -       .       ID=chain_31;Target=asmbl_31 36 94 +
#373     alignAssembly-ME49_genes_tgrh88_pasa    cDNA_match      161728  161762  .       -       .       ID=chain_31;Target=asmbl_31 1 35 +

## OUTPUT

#tgcast.assembly.1	pasa	intron	23970	24753	.	+	.	chain_1_1;source=M
#tgcast.assembly.1	pasa	intron	24811	25743	.	+	.	chain_1_2;source=M
#tgcast.assembly.1	pasa	intron	25900	26414	.	+	.	chain_1_3;source=M
#tgcast.assembly.1	pasa	intron	26590	27190	.	+	.	chain_1_4;source=M
#tgcast.assembly.1	pasa	intron	27297	27708	.	+	.	chain_1_5;source=M

########################################

my $usage = "$0 -i <pasa_gff3_file> -d database_name\n\n";
my %arg = @ARGV;
die $usage unless $arg{-i} and $arg{-d};
open (GFF3, "$arg{-i}") || die "ERROR, I cannot open $arg{-i}: $!\n\n";
my (%counter, $chain_id, $old_chain_id, $prev_line);
while(<GFF3>){
	chomp;
	next unless m/cDNA_match/;
	my ($asmbl, undef,undef, $end5, $end3, undef, $strand, undef, $comment) = split /\t/;
	my $chain_id = $1 if $comment =~ m/ID=chain_(\d+)/;
	die "ERROR, I cannot fetch chain_ID from $comment\n\n" unless $chain_id;
	$end3++;
	if($chain_id != $old_chain_id){
		## new transcript
		$prev_line = "$arg{-d}.assembly.$asmbl\tpasa\tintron\t$end3";
	}
	else {
		$counter{$chain_id}++;
		my $com = "chain_".$chain_id."_".$counter{$chain_id}.";source=M";
		$end5--; 
		print "$prev_line\t$end5\t.\t$strand\t.\t$com\n";
		$prev_line = "$arg{-d}.assembly.$asmbl\tpasa\tintron\t$end3";
			
	}
	$old_chain_id = $chain_id;
}
close GFF3;



