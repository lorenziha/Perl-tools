#!/usr/local/bin/perl
use strict;
use Data::Dumper;

my $usage = "$0 -a <DAGchainer_aligncoords_file>\n";

die $usage if !@ARGV || $ARGV[0] eq '-h';
my %arg = @ARGV;
die $usage unless ($arg{-a});

open (ALIGN, "<$arg{-a}") || die "ERROR, I cannot open $arg{-a}: $!\n\n";
my @file = <ALIGN>;
my @file_copy = @file;
close ALIGN;

my (%query, %subject);
while( @file){
	my $syntenic_block = &get_next_syntenic_block(@file);
	my $query = &get_features($syntenic_block, 'query');
	my $subject = &get_features($syntenic_block, 'subject');
	push @{ $query{ $query->{asmbl} } }, $query;
	push @{ $subject{ $subject->{asmbl} } }, $subject;
}

## Sort %query;
my %query_ID;
foreach my $key (keys %query){
	my @unsorted = @{$query{$key}};
	my @sorted = sort { $a->{mid} <=> $b->{mid}} @unsorted;
	@sorted = &set_ends(@sorted);
	foreach my $entry (@sorted){
		$query_ID{ $entry->{id} } = $entry;
	}
	@{$query{$key}} = @sorted;
	
}

## Sort %subject
my %subject_ID;
foreach my $key (keys %subject){
	my @unsorted = @{$subject{$key}};
	my @sorted = sort { $a->{mid} <=> $b->{mid}} @unsorted;
	@sorted = &set_ends(@sorted);
	foreach my $entry (@sorted){
		$subject_ID{ $entry->{id} } = $entry;
	}
	@{$subject{$key}} = @sorted;
}

@file = @file_copy;
while( @file){
	my $syntenic_block = &get_next_syntenic_block(@file);
	my $query = &get_features($syntenic_block, 'query');
	my $subject = &get_features($syntenic_block, 'subject');
	my ($q_start, $q_end, $s_start, $s_end) = ('CONT','CONT','CONT','CONT');
	if ( $query_ID{ $query->{id} }->{block_start} == 1){$q_start = 'START'}
	if ( $query_ID{ $query->{id} }->{block_end} == 1){$q_end = 'END'}
	if ( $subject_ID{ $subject->{id} }->{block_start} == 1){$s_start = 'START'}
	if ( $subject_ID{ $subject->{id} }->{block_end} == 1){$s_end = 'END'}
	
	if ($syntenic_block->[0] =~ m/reverse/){
		## invert $q_start, $q_end values
		($q_start,$q_end) = ($q_end,$q_start);
	}
	
	#print Dumper($syntenic_block);
	
	## print header
	my $header = shift @{$syntenic_block};
	print "$header";
	
	## print block start status
	print "$s_start($subject->{asmbl})-$q_start($query->{asmbl})\n";
	
	## print remaining of the block
	print @{$syntenic_block};
	
	## print block end status
	print "$s_end($subject->{asmbl})-$q_end($query->{asmbl})\n";
}

exit(0);

###########################################s###
sub extract_ID {
	my @block = @_;
	my %ID_hash;
	foreach my $entry (@block){
		$ID_hash{ $entry->{id} } = $entry;
	}
	return \%ID_hash;
}

sub set_ends {
	my @block = @_;
	$block[0]->{block_start} = 1;
	$block[-1]->{block_end} = 1;
	return @block;
}

sub get_features {
	## get query start, end, ID 
	my $block = shift;
	my $type = shift;
	my $start_line = $block->[1];
	my $end_line = $block->[-1];
	my ($start,$end,$mid_point,$asmbl,$id) = &process_start_end_lines($type,$start_line,$end_line);
	return { 	start => $start,
				end => $end,
				mid => $mid_point,
				asmbl => $asmbl,
				id => $id,
				block_start => 0,
				block_end => 0,
			};
}

sub process_start_end_lines {
	my ($type, $start_line,$end_line) = @_;
	my $padding = 0;
	if ($type eq 'query'){
		$padding = 4;
	}
	my @start_line = split(m/[\s\t]+/,$start_line);
	my @end_line = split(m/[\s\t]+/,$end_line);
	my $mid_point = int( ($start_line[$padding + 2] + $end_line[$padding + 2]) / 2 );
	my $id = $start_line[$padding].'_'.$start_line[$padding + 2].'_'.$end_line[$padding + 2];
	return ($start_line[$padding + 2], $end_line[$padding + 2], $mid_point, $start_line[$padding], $id);
}

sub get_next_syntenic_block {
	my @block;
	push @block, shift @file;
	while ($file[0] =~ m/^\d/){
		push @block, shift @file;
	}
	return \@block;
}
