#!/usr/local/bin/perl -w
use strict;
use Bio::Tools::Run::Phylo::PAML::Codeml;
use Bio::Tools::Run::Alignment::Clustalw;
use Data::Dumper;
 
# for projecting alignments from protein to R/DNA space
use Bio::Align::Utilities qw(aa_to_dna_aln);
# for input of the sequence data
use Bio::SeqIO;
use Bio::AlignIO;


my $usage = "$0 -f fasta -o outputdir\n\n";
my %arg = @ARGV;
die $usage unless $arg{-f} && $arg{-o};

my $aln_factory = Bio::Tools::Run::Alignment::Clustalw->new('-quiet' => 1);
my $seqdata = $arg{-f};
my $dir = $arg{-o};
unless (-d $dir){
	die "ERROR, I cannot find $dir\n\n";
}
 
my $seqio = new Bio::SeqIO(-file   => "$dir/$seqdata",
                           -format => 'fasta');
my %seqs;
my @prots;
# process each sequence
while ( my $seq = $seqio->next_seq ) {
    $seqs{$seq->display_id} = $seq;
    # translate them into protein
    my $protein = $seq->translate();
    my $pseq = $protein->seq();
    if( $pseq =~ /\*/ &&
        $pseq !~ /\*$/ ) {
          warn("provided a CDS sequence with a stop codon, PAML will choke!");
          exit(0);
    }
    # Tcoffee can't handle '*' even if it is trailing
    $pseq =~ s/\*//g;
    $protein->seq($pseq);
    push @prots, $protein;
}
 
if( @prots < 2 ) {
    warn("Need at least 2 CDS sequences to proceed");
    exit(0);
}
 
#open(OUT, ">align_output.txt") ||  die("cannot open output align_output for writing");
# Align the sequences with clustalw
my $aa_aln = $aln_factory->align(\@prots);

# get tree
my $tree = $aln_factory->tree(\@prots);

# project the protein alignment back to CDS coordinates
my $dna_aln = aa_to_dna_aln($aa_aln, \%seqs);
 
my @each = $dna_aln->each_seq();
 
my $kaks_factory = Bio::Tools::Run::Phylo::PAML::Codeml->new
                   ( -params => { 'runmode' => 0,
                                  'seqtype' => 1,
				  'model' => 0,
				  ,	
                                } );
 
# cleanup results ? 0=yes; 1=no
$kaks_factory->save_tempfiles(1);

# set the alignment object
$kaks_factory->alignment($dna_aln);

# set the tree object
$kaks_factory->tree($tree);
#print "kaks_factory=".Dumper($kaks_factory);
# run the KaKs analysis
my ($rc,$parser) = $kaks_factory->run();

#exit(0);

my $result = $parser->next_result;
my $MLmatrix = $result->get_MLmatrix();

#get tmp working directory
my $tmpdir = $kaks_factory->tempdir();  
my @otus = $result->get_seqs();

open(MLC,"$tmpdir/mlc") || die "ERROR, I cannot open $tmpdir/mlc file: $!\n\n";
while(<MLC>){
	chomp;
	if(m/^omega\s+\(dN\/dS\)\s=\s+(\d+\.\d+)/){
		print "$seqdata\tdN/dS\t$1\n";
		last;
	}
}
close MLC;

`cp $tmpdir/mlc $dir/$seqdata.mlc`;
`rm -Rf $tmpdir`;




# this gives us a mapping from the PAML order of sequences back to
# the input order (since names get truncated)
#my @pos = map {
#    my $c= 1;
#    foreach my $s ( @each ) {
#      last if( $s->display_id eq $_->display_id );
#      $c++;
#    }
#    $c;
#} @otus;


#print "Result= ".Dumper($result);

#print "\n\notus=".Dumper(@otus);

#print OUT join("\t", qw(SEQ1 SEQ2 Ka Ks Ka/Ks PROT_PERCENTID CDNA_PERCENTID)),"\n";
#foreach my $i ( 0 .. $#otus -1 ) {
#  foreach my $j ( $i+1 .. $#otus ) {
#    my $sub_aa_aln  = $aa_aln->select_noncont($pos[$i],$pos[$j]);
#    my $sub_dna_aln = $dna_aln->select_noncont($pos[$i],$pos[$j]);
#    print OUT join("\t", $otus[$i]->display_id,
#                         $otus[$j]->display_id,$MLmatrix->[$i]->[$j]->{'dN'},
#                         $MLmatrix->[$i]->[$j]->{'dS'},
#                         $MLmatrix->[$i]->[$j]->{'omega'},
#                         sprintf("%.2f",$sub_aa_aln->percentage_identity),
#                         sprintf("%.2f",$sub_dna_aln->percentage_identity),
#                         ), "\n";
#  }
#}
