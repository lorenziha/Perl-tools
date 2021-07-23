#!/usr/local/bin/perl -w
use strict;

my $usage = "$0 -w window_size -r RepeatMasker_out_file -n T/F collapse_same_repeats [F]\n\n";

die $usage if scalar @ARGV == 0 or $ARGV[0] eq '-h';
my %arg = @ARGV;

my $file = $arg{-r};
my $window = $arg{-w};
my $NO_DUPLICATES = $arg{-n} eq 'T' ? 1 : 0;

my %rm; # hash with RepeatMasker output
my %repeat_types; # hash containing as keys all the different repeat types

## Load repeats
open(FHI,"<$file") || die "ERROR, I cannot open $file: $!\n";
while(<FHI>){
	chomp;
	next if m/^[\s\t]*$/;
	s/^\s+//;

	my $repeat = &get_repeat_info( $_ );
	push @{ $rm{ $repeat->{asmbl} } }, $repeat;
	$repeat_types{ $repeat->{repeat} }++;
}
close FHI;

my %cluster_counter;

## Look for cluster types
foreach my $asmbl ( keys %rm ){
	my @cluster_list = &get_clusters( @{ $rm{$asmbl} } );
	
	## Count clusters
	foreach my $cluster_id (@cluster_list){
		$cluster_counter{$cluster_id}++;
	}
}

## Print results
print "INPUT_FILE= $file\n";
print "WINDOW= $window\n\n";
my %sizes;
foreach my $cluster_id (sort {length($a) <=> length($b)} keys %cluster_counter){
	my @current_cluster = split('-', $cluster_id);
	my $number_of_members = scalar( @current_cluster );
	$sizes{ $number_of_members } += $cluster_counter{ $cluster_id };
	print "Cluster_ID= $cluster_id\t$cluster_counter{$cluster_id}\n" if $cluster_id;
}

print "\n\nCluster_size\tfrequency\n";

foreach my $cluster_size (sort {$a <=> $b} keys %sizes ){
	next unless $cluster_size;
	print "$cluster_size\t$sizes{$cluster_size}\n";
}

## Count number of times a pair of two repeats cluster together
my @repeat_types = keys %repeat_types;
my %repeat_pair_counter; # hash containing counts of all possible pair of repeats

print "\nRepeat_pair\tcounts\n";

while ( my $repeat1 = shift @repeat_types){
	foreach my $repeat2 (@repeat_types){
		my $repeat_pair = $repeat1.'-'.$repeat2;
		$repeat_pair_counter{ $repeat_pair } = 0; #initialize hash for each pair
		
		foreach my $cluster_id ( keys %cluster_counter ){
			## check if cluster_id contains both repeats
			if ($cluster_id =~ m/$repeat1/ && $cluster_id =~ m/$repeat2/ ){
				$repeat_pair_counter{ $repeat_pair } += $cluster_counter{$cluster_id};
			}
		}
		
		print "$repeat_pair\t$repeat_pair_counter{ $repeat_pair }\n";	
	}
}

#############################################################
sub get_clusters {
	my @asmbl = @_;
	my $old_hit = shift @asmbl;
	my @block = ($old_hit->{repeat});
	my @clusters; #stores list of blocks
	
	## if just a single hit in this asmbly
	push @clusters, $old_hit->{repeat} unless scalar @asmbl;
	if(scalar @asmbl){
	foreach my $new_hit (@asmbl){
			if( ($new_hit->{end5} - $old_hit->{end3}) <= $window ){
				
				## ADD ELEMENT TO CURRENT BLOCK
				
				push @block, $new_hit->{repeat};
			}
			if ( ( ($new_hit->{end5} - $old_hit->{end3}) > $window) 
					 or $new_hit eq $asmbl[-1] ){
				
				## IF DISTANCE > WINDOW OR LAST ELEMENT
				## SAVE OLD BLOCK AND START A NEW BLOCK
				
				@block = sort {$a cmp $b} @block;
				
				## transform cluster A-B-A into A-B
				@block = &remove_duplicates(@block) if $NO_DUPLICATES;
				my $cluster_id = join(' /-/ ',@block);
				@block = (); # reset @block
				push @clusters, $cluster_id;
				
				## Add last element of the asmbl as new cluster if it doesn't cluster w/ old_hit
				push @clusters, $new_hit->{repeat} if ( ($new_hit->{end5} - $old_hit->{end3}) > $window);
			}

			$old_hit = $new_hit;
		}
	}
	return @clusters;
}

sub remove_duplicates {
	my @redundant = @_;
	my %non_redundant;
	@non_redundant{@redundant} = @redundant;
	return keys %non_redundant;
}

sub get_repeat_info {
	my $line = shift;
	my @row = split (m/[\s\t]+/,$line);
	return {	asmbl => $row[4],
				end5  => $row[5],
				end3  => $row[6],
				strand => $row[8],
				repeat => $row[9]
				}; 
}
