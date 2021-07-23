#!/usr/local/bin/perl
use strict;
use Bio::SeqIO;
use Bio::Seq;

my $usage = "$0 -i <pfma.htab file> -f <pep_fasta_file>\n\n";
my %arg = @ARGV;
die $usage unless $arg{-i} or $arg{-f};

## get protein length info fro fasta file
my $inseq = Bio::SeqIO->new(
                            -file   => "<$arg{-f}",
                            -format => "fasta",
                            );
my %prot_len;
while (my $seq = $inseq->next_seq) {
        my $pub_locus = $seq->primary_id;
        my $seqlen  = $seq->length;
        $prot_len{$pub_locus} = $seqlen;
	#print "$pub_locus = $seqlen\n";
}

## get hmm hit values
open (HTAB,"$arg{-i}") || die "ERROR, I cannot open $arg{-i} : $!\n\n";
my %protein;
my %proteins_with_pfams;
while(<HTAB>){
	my @x = split /\t/;
	my ($pfam, $plen, $pub_locus,$pstart, $pend, $start, $end, $score, $tc, $evalue_dom, $description ) 
	= ($x[0],$x[2], $x[5], $x[6], $x[7], $x[8], $x[9], $x[11], $x[17],$x[20],$x[15]); 
	next if $score < $tc;
	my $pfam_coverage = ($pend + 1 - $pstart) / $plen;
	## next unless $pfam_coverage >= 0.5; ## discard pfam hits shorter than 60% of the pfam domain
	$pub_locus =~ s/\s+//g;
	my $teval = $evalue_dom == 0 ? (1/1e-100) : (1/$evalue_dom);
	#print "teval = $teval\n";
	my $prot_len = $prot_len{$pub_locus} || warn "WARNING, no length info for -$pub_locus-\n\n";
	push @{ $protein{$pub_locus} }, [ $pfam, $plen, $prot_len, $pub_locus, $start, $end, $teval, $description, $evalue_dom  ];
	$proteins_with_pfams{$pub_locus}++;
}
close HTAB;

print "## PFAM domains have not to overlap with other domain with lower e-value and have a total score above the truscted cutoff\n";
foreach my $protein ( keys %proteins_with_pfams ){
	my @array = @{ $protein{$protein} };
	my @sorted = sort { $b->[6] <=> $a->[6]  } @array;
	my $seq = 'P' x $sorted[0]->[2];
	my $overlap;
	foreach my $hits (@sorted){
		my ($pfam, $plen, $prot_len, $pub_locus, $start, $end, $teval, $description, $evalue) = @{$hits};
		my $hit_span = 'N' x (($end - $start) + 1);
		my $hit = substr($seq, $start, ($end - $start) + 1);
		if ($hit =~ m/N/){
			## overlaps another domain
			$overlap++;
		}
		else {
			substr($seq, $start, ($end - $start) + 1, $hit_span);
			print "$pub_locus\t$pfam\t$evalue\t$start - $end\t$description\n";
		}
	}
}




