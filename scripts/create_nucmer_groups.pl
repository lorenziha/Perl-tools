#!/usr/local/bin/perl
use strict;
use lib '/home/hlorenzi/PERL/SCRIPTS';
use compare_hits;

my $usage ="create_nucmer_groups.pl -i <nucmer_coord_file> -t <#, minimum hit size to consider> [-max <#, maximum hit size [1Mb]> -s minimum_asmbl_size [0] -id minimum % of identity [0] -c coverage between hits [50]]\n";

die $usage if $ARGV[0] eq '-h';
my %arg = @ARGV;

my %asmbl_list;
my $group_id;
my %asmbl_group;
my $asmbl_coords_ref;
my %coord_5;
my %coord_3;
my $len = 7000;
my $asmbl_pair;
my $min_size = $arg{-t} || die $usage;
my $max_size = $arg{-max} || 1000000;
my $min_asmbl_len = $arg{-s} || 0;
my %column_number;
my $id_cutoff = $arg{-id} || 0;
my $coverage = $arg{-c} || 50;

open(FHI,"<$arg{-i}") || die "I cannot find $arg{-i}:$!\n";

ASMBL_PAIR:
while(<FHI>){
	
	chomp;
	next ASMBL_PAIR if m/^[\s\t]*$/;
	next if scalar( split m/\s+/) < 7;
	s/[\|]//g;
	
	# get column positions
	if(m/\[S1\]/){
		
		s/\[TAGS\]/\[REF\]   [QRY\]/; 
		my $column_count = 0;
		foreach my $label (split m/\s{2,}|\t/){
			$column_number{$label} = $column_count;
			$column_count++;
		}
		die "Coords file must contain asmbl length info!\n" if (!exists $column_number{'[LEN R]'});
		die "Coords file must contain % ID info!\n"         if (!exists $column_number{'[% IDY]'});
	}
	next unless m/^\s+\d+\s+\d+\s+\d+\s+\d+/;
	#get coords and asmbls

	my @row = split m/[\s\t]+/;
	my ($ref_start, $ref_end, $q_start, $q_end, $hit_len, $perc_id, $ref_len, $q_len, $reference, $query) 
		= ($row[1],$row[2],$row[3],$row[4],$row[ $column_number{'[LEN 1]'} ], $row[ $column_number{'[% IDY]'} ], $row[ $column_number{'[LEN R]'} ],$row[ $column_number{'[LEN Q]'} ],$row[-2],$row[-1]);
	
	#print "processing... $ref_start, $ref_end, $q_start, $q_end, $hit_len, $perc_id, $ref_len, $q_len, $reference, $query\n";
	if( $reference =~ m/\w+\d\.\w+\.(\d+)\.\d+/){ $reference = $1;}
	if( $query =~ m/\w+\d\.\w+\.(\d+)\.\d+/){ $query = $1;}
	
	next ASMBL_PAIR unless ( $hit_len >= $min_size );
	next ASMBL_PAIR unless ( $hit_len <= $max_size );
	next ASMBL_PAIR unless ( $ref_len >= $min_asmbl_len and $q_len >= $min_asmbl_len );
	next ASMBL_PAIR unless ( $perc_id >= $id_cutoff);
	next ASMBL_PAIR if ( $reference == $query); ## avoid self hit
	next ASMBL_PAIR unless ( $reference =~ m/^\d+$/ and $query =~ m/^\d+$/);
	
	print STDERR "processing... ref5=$ref_start, ref3=$ref_end, q5=$q_start, q3=$q_end, hit_len=$hit_len, \%id=$perc_id, ref_len=$ref_len, q_len=$q_len, ref_id=$reference, q_id=$query\n";

	
	# load data into asmbl_pair hash
	$asmbl_pair->{$reference}{ref_len}   = $ref_len;
	$asmbl_pair->{$query}{q_len}         = $q_len;

	
	unless ($ref_start > 0 and $ref_end > 0 and $q_start > 0 and $q_end > 0 ){
		print STDERR "Coord error!! : ref_5 $ref_start ref_3 $ref_end q_5 $q_start q_3 $q_end\n";
		sleep 2;
	}
	
	# sort hit coordinates
	( $coord_5{$reference},$coord_3{$reference} ) =  ($ref_start,$ref_end);
	( $coord_5{$query}    ,$coord_3{$query}     ) =  ($q_start  ,$q_end  );
	
	# if both asmbls already loaded
	if(exists $asmbl_list{$reference} and exists $asmbl_list{$query}){
		
		# check if both asmbls belong to the same group
		my @groups_to_fuse = ();
		my $number_of_hits = 0;
		
		GROUP:
		foreach my $group_id ( sort{ $asmbl_group{$a} <=> $asmbl_group{$b} } keys %asmbl_group){
			$number_of_hits = scalar( grep( m/^($reference|$query)$/, @{ $asmbl_group{$group_id} } ) );
			my @hits = grep( m/^$reference|$query$/, @{ $asmbl_group{$group_id} } );
			# asmbls not present in current group
			next GROUP if $number_of_hits == 0;
			
			# both asmbl in the group, do nothing and take another pair
			if ($number_of_hits == 2){
				
				# compare reference coords with reference group coords to see if they overlap
				my ($group_ref_5,$group_ref_3) = split(m/,/,$asmbl_coords_ref->{$reference}{$group_id});
				#print "A1= ref=$reference, group_id=$group_id, $coord_5{$reference},$coord_3{$reference},$group_ref_5,$group_ref_3\n";
				my $coord_analysis = compare_hits($coord_5{$reference},$coord_3{$reference},$group_ref_5,$group_ref_3);
				
				if ( $coord_analysis->{coverage1} >= $coverage or $coord_analysis->{coverage2} >= $coverage ){
					next ASMBL_PAIR;
				}
				
				# if the hit coordinates do not overlap with at least 50% of coverage
				else {
					$number_of_hits = 0;
				}
			}			
			
			# each asmbl in a different group => fuse groups
			elsif($number_of_hits == 1){
				my ($asmbl_already_present) = grep( m/^($reference|$query)$/, @{ $asmbl_group{$group_id} } );
				
				# compare reference coords with reference group coords to see if they overlap
				my ($group_ref_5,$group_ref_3) = split(m/,/,$asmbl_coords_ref->{$asmbl_already_present}{$group_id});
				#print "A2= asmbl=$asmbl_already_present, group_id=$group_id, $coord_5{$reference},$coord_3{$reference},$group_ref_5,$group_ref_3\n";
				my $coord_analysis = compare_hits($coord_5{$asmbl_already_present},$coord_3{$asmbl_already_present},$group_ref_5,$group_ref_3);
				
				if ( $coord_analysis->{coverage1} >= $coverage or $coord_analysis->{coverage2} >= $coverage ){
					
					push @groups_to_fuse, $group_id;
				}

			}
			elsif($number_of_hits > 2){
				print join( ',' , @hits ),"\n";
				die "Exception [4]: more than 2 hits per group ($number_of_hits)\n";
			}			
		}
		
		# check how many groups to fuse
		my $number_of_groups = scalar(@groups_to_fuse);
		die "Exception [2]: number of goups to fuse is more than 2 ($number_of_groups)\n" if $number_of_groups > 2;
	
		if ($number_of_groups == 2){
		
			# fuse groups
			@groups_to_fuse = sort @groups_to_fuse;
			push @{ $asmbl_group{ $groups_to_fuse[0] } }, @{ $asmbl_group{ $groups_to_fuse[1] } };
			@{ $asmbl_group{ $groups_to_fuse[0] } } = remove_duplicated_entries( @{ $asmbl_group{ $groups_to_fuse[0] } } );
		
			# assign old coords hash ref to new group
			foreach my $asmbl ( @{ $asmbl_group{ $groups_to_fuse[1] } } ){
				$asmbl_coords_ref->{$asmbl}{$groups_to_fuse[0]} = $asmbl_coords_ref->{$asmbl}{$groups_to_fuse[1]};
			}
		
			# remove old group 
			delete $asmbl_group{ $groups_to_fuse[1] };
		}
		elsif ($number_of_groups == 1){
			# remove $query from the asmbl list to make next condition true (see next if loop)
			delete $asmbl_list{$query};
		}
		
	}	
	
	# Add missing asmbl to a group
	if( (exists $asmbl_list{$reference} and !exists $asmbl_list{$query} )
		or (!exists $asmbl_list{$reference} and exists $asmbl_list{$query} )
			 ){
		
		GROUP_B:
		foreach my $group_id ( sort{ $asmbl_group{$a}<=>$asmbl_group{$b} } keys %asmbl_group){
			
			if( my ($present_asmbl) = grep( m/^($reference|$query)$/, @{ $asmbl_group{$group_id} } ) ){
				my $new_asmbl;
		
				if($reference == $present_asmbl){
					$new_asmbl = $query;
				}
				elsif($query == $present_asmbl){
					$new_asmbl = $reference;
				}
				else {
					die "Exception [1]: unknown asmbl_id ($present_asmbl)\n";
				}
				
				# compare reference coords with reference group coords to see if they overlap
				my ($group_ref_5,$group_ref_3) = split(m/,/,$asmbl_coords_ref->{$present_asmbl}{$group_id});
				#print "B= asmbl=$present_asmbl, group_id=$group_id, $coord_5{$reference},$coord_3{$reference},$group_ref_5,$group_ref_3\n";
				my $coord_analysis = compare_hits($coord_5{$present_asmbl},$coord_3{$present_asmbl},$group_ref_5,$group_ref_3);
				
				if ( $coord_analysis->{coverage1} >= 50 or $coord_analysis->{coverage2} >= 50 ){
					push @{ $asmbl_group{$group_id} }, $new_asmbl;
					@{ $asmbl_group{ $group_id } } = remove_duplicated_entries( @{ $asmbl_group{ $group_id } } );
					$asmbl_list{$new_asmbl} = 1;
					$asmbl_list{$present_asmbl} = 1;
					
					# create new coord hash ref for new asmbl
					$asmbl_coords_ref->{$new_asmbl}{$group_id} = "$coord_5{$new_asmbl},$coord_3{$new_asmbl}";
					
					last GROUP_B;
				}
				else {
					# remove asmbl_ids from the list of asmbls
					delete $asmbl_list{$new_asmbl};
					delete $asmbl_list{$present_asmbl};
					next GROUP_B; # No hit overlap => do something else, create new group?
				}
			}
		}
	}
	
	# Create new group with both asmbls
	elsif(!exists $asmbl_list{$reference} and !exists $asmbl_list{$query}){
		# incorporate asmbl ids to the list of asmbls
		( $asmbl_list{$reference}, $asmbl_list{$query} ) = (1,1);
		
		# create new group
		$group_id++;
		push @{ $asmbl_group{$group_id} }, ($reference,$query);
		
		# create new coord hash ref
		$asmbl_coords_ref->{$reference}{$group_id} = "$coord_5{$reference},$coord_3{$reference}";
		$asmbl_coords_ref->{$query    }{$group_id} = "$coord_5{$query},$coord_3{$query}";
		
	}
}

# PRINT RESULTS
my %asmbl_pair;

foreach my $group_id ( sort{ $a<=>$b } keys %asmbl_group){
	my @current_group = ();
	my $number_of_members = 0;
	foreach my $asmbl ( @{ $asmbl_group{$group_id} } ){
		$number_of_members++;
		my $asmbl_coord = $asmbl."[$asmbl_coords_ref->{$asmbl}{$group_id}]";
		push @current_group,$asmbl_coord;
		$asmbl_pair{"$asmbl_group{$group_id}->[0]".'_'."$asmbl"} = 1;
	}
	
	print 'Group '.$group_id."\t$number_of_members\t".join(':', @current_group),"\n";
}

close FHI;

my $count = 0;
my $total_len_wo_1_copy = 0;
my $total_len = 0;

foreach my $group_id ( sort{ $a<=>$b } keys %asmbl_group){
	my @current_group = ( @{ $asmbl_group{$group_id} } );
	my $pair_id;
	$count++;
	print "\n$count) Group $group_id\t".join(',', @current_group),"\n";
	foreach my $asmbl ( @current_group ){

		my $strand;
		
		# if first element of the array => reference seq
		if ($asmbl == $current_group[0]){
		
			$pair_id =  $current_group[0]."_".$current_group[1];



			my ($ref_start,$ref_end) = split(/,/,$asmbl_coords_ref->{$asmbl}{$group_id});
			my $ref_len = $asmbl_pair->{$asmbl}{ref_len};
			
			
			if( ($ref_end - $ref_start ) > 0){
				$strand = "[1]->[$ref_start,$ref_end]->[$ref_len]";
				
			}
			else {
				$strand = "[$ref_len]<-[$ref_start,$ref_end]<-[1]";
				
			}

		}
		
		# if other array element => query sequence
		else {
			$pair_id =  $current_group[0]."_".$asmbl;
			my ($q_start,$q_end) = split(/,/,$asmbl_coords_ref->{$asmbl}{$group_id});
			my $q_len = $asmbl_pair->{$asmbl}{q_len};
			
			
			if( ($q_end - $q_start) > 0){
				$strand = '[1]->['.$q_start.','.$q_end."]->[$q_len]";
				$total_len_wo_1_copy += $q_end - $q_start;
				$total_len += $q_end - $q_start;
			}
			else {
				$strand = '['.$q_len.']<-['.$q_start.','.$q_end."]<-[1]";
				$total_len_wo_1_copy += $q_start - $q_end;
				$total_len += $q_start - $q_end;
			}
		}
		print "$asmbl:\t$strand\n";
	}
	

}


sub remove_duplicated_entries {
	my %array_wo_duplicates;
	foreach my $entry (@_){
		$array_wo_duplicates{$entry} = $entry;
	}
	return keys %array_wo_duplicates;
}

print "\n***********************************\nTotal length without 1 copy: $total_len_wo_1_copy\nTotal length: $total_len\n***********************************\n";
