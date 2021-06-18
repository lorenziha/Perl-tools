#!/Users/lorenziha/bin/miniconda3/bin/perl
use strict;

my $usage = "$0 -i <fasta file> -w <window_size [120bp]>\n\n";
my %arg = @ARGV;
my $WINDOW = $arg{-w} || 120;
die $usage unless $arg{-i};

open (FASTA, "<$arg{-i}") || die "ERROR, I cannot open $arg{-i}: $!\n\n";
my ($id, $seq, %seq, @ids);
while(<FASTA>){
	chomp;
	next if m/^[\s\t]*$/;
	if(m/^>(\S+)/){
		#header
		$id = $1;
		push @ids, $id;
	}
	elsif(m/^([acgtn]+)$/i){
		#seq
		$seq{$id} .= $1;
	}
}
close FASTA;

# Calculate GC content
# Output format has to be chr_id<tab>position<tab>gc_content
foreach my $id (@ids){
	my $seq = $seq{$id};
	for (my $pos = 0; $pos < length($seq)-$WINDOW; $pos++) {
		my $subseq = substr($seq,$pos,$WINDOW);
		my $adjusted_pos = $pos + int($WINDOW/2);
		my $GC = &get_gc($subseq);
		print "$id\t$adjusted_pos\t$GC\n";
	}
}
exit(0);

##############################################################
sub get_gc{
	my $subseq = shift;
	my $total_len = length($subseq);
	$subseq =~ s/[AT]//gi;
	my $gc_len = length($subseq);
	my $GC_perc = int($gc_len * 100 / $total_len);
	return ($GC_perc);
}