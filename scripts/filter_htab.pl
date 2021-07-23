#!/usr/local/bin/perl 
=head1 NAME

filter_htab.pl - Filters htab files by e-value, coverage, and trusted/gathering cutoffs

=head1 USAGE

	filter_htab.pl -i XXXX.htab file
	
=head1 OPTIONS

    REQUIRED:

    -i input htab list file

	OPTIONS
	
	-e e-value cutoff for whole match [default 0.001]
	-c coverage cutoff 0-100 [default 100]
	-t use trusted cutoff for whole score match**
	-g use gathering cutoff for whole score match**
	-h prints this help
	
	** if both, then uses trusted cutoff

=head1  OUTPUT

The program prints to the following output file:

	input_file.htab_eXXcovXXcut[tg-].htab
	eXXcovXXcut[tg-].htab.list
		
=head1  CONTACT

    Hernan Lorenzi
    hlorenzi@jcvi.org or
    
=begin comment

    ## legal values for status are active, inactive, hidden, unstable
    status: active
    keywords: molecular weight isoelectric pI pepstats

=end comment

=cut

use Pod::Usage;
use Getopt::Std;
use strict;
use Data::Dumper;

my $opt = {};
&getopts('hi:e:c:tg', $opt);

if($opt->{h}) { &pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} ); }

my $htab_list_file = $opt->{i} || &pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
my $evalue = $opt->{e} || 0.001;
my $coverage = $opt->{c} || 100;
if ( $coverage > 100 || $coverage < 0 ) {&pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );}
unless (-e $htab_list_file) { die "ERROR, I cannot find $htab_list_file"; }
my $cutoff_type = $opt->{t}? 'T' : $opt->{g}? 'G' : '-';
my $tag = 'e_'.$evalue.'_cov_'.$coverage.'_cut_'.$cutoff_type.'.htab';
my $list_file = $tag.'.list';

open (LIST,"<$htab_list_file") || die "ERROR, I cannot open $htab_list_file: $!\n\n";
while(<LIST>){
	chomp;
	next if m/^[\s\t]*$/;
	my $htab_file = $_;
	warn "WARNING, I cannot find $htab_file !!\n\n", next if (! -e $htab_file);
	my $output_file = $htab_file."_$tag";

	
	## read htab file
	open (HTAB,"<$htab_file") || die "ERROR, I cannot open $htab_file: $!\n\n";
	open (OFH, ">$output_file")|| die "ERROR, I cannot open $output_file: $!\n\n";
	while(<HTAB>) {
		chomp;
		my $line = $_;
		my $info = &get_htab_info($line,$cutoff_type);
		#print Dumper($info);

		if ($info->{evalue} <= $evalue && $info->{coverage} >= $coverage && $info->{score} >= $info->{cutoff} ) {

			print  OFH "$line\n";
		}
	}
	close HTAB;
	close OFH;

	## add file to list file
	open (OUTLIST, ">>$list_file") || "ERROR, I cannot open $list_file: $!\n\n";
	print OUTLIST "$output_file\n";
	close OUTLIST;
}
close LIST;
exit(0);

##############################################################################
sub get_htab_info {
	my ($line,$cutoff_type) = @_;
	my @line = split(m/\t/,$line);
	my ($end5, $end3, $hmm_len) = ($line[6], $line[7], $line[2]);
	$hmm_len = 1 unless $hmm_len > 0; ## fixes the problem that some domain lengths are missing from the htab files
	my $coverage = ( abs($end3-$end5) + 1) * 100 / $hmm_len;
	my $cutoff = $cutoff_type eq 'T' ? $line[21] : $cutoff_type eq 'G' ? $line[23] : -1000000;
	return {	evalue => $line[19],
				coverage => $coverage,
				cutoff => $cutoff,
				score => $line[12],
				};
}
