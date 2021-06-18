#!/usr/bin/perl
#
my $usage = "$0 -k <keys_file> [-ki <key file index [0]> ; -d <data_file> | STDIN ; -di <data index [0]> ;-m M|MM [match|mismatch [MM]] ; -h help]\n\n";

my %arg = @ARGV;
die $usage if $arg{-h} || !exists($arg{-k});
my $RINDEX = $arg{-ki} || 0;
my $QINDEX = $arg{-di} || 0;
my $MISMATCH = $arg{-m} eq "M"? 0 : 1;

# Load reference ids to query in the input file
print STDERR "Loading keys...\n";
open (REF, "<$arg{-k}") || die "I cannot find file $arg{-d}: $!\n\n";
my %h;
while(<REF>){
	chomp;
	my @x=split /\t/;
	$h{ $x[$RINDEX] }++;
}
close REF;

print STDERR "Reading data file...\n";
my $pace = 1000000; #controls progress bar
my ($QUERY);

# Do the query
if ($arg{-d}){
	open( $QUERY , "$arg{-d}") || die "I cannot find file $arg{-d}: $!\n\n";
} else {
	$QUERY = *STDIN; 
}
my $counter = 0;
while(<$QUERY>){
	chomp;
	$counter++;
	my @x = split /\t/;
	if( $MISMATCH ){
		print "$_\n" unless $h{ $x[$QINDEX] };	
	} else {
                print "$_\n" if $h{ $x[$QINDEX] };
	}
	if ($counter >= $pace){
		print STDERR "=";
		$counter = 0;
	}
}
print STDERR "\n";
close $QUERY if $arg{-d};
