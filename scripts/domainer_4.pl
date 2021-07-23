#!/usr/local/bin/perl

use DBI;
use strict;
use CGI qw/:standard/;


$|=1;

      
# hlorenzi Oct 14, 2005
# domainer_4

## Variable definitions
my $version = 'v4.0';
my $dbh;
my @filedom;
my $e_val = 0.001;
my $total_score = "gathering_cutoff";
my %db;
my $accession;
my @filein = ();
my @filedb =();
my $feat_name;
my $hmm_com_name;
my $is_public = 1;

## Parameters entry

my %arg = @ARGV;

## INPUT FILES/DATABASES

if( exists( $arg{-h} ) ){ die Help() }
if( exists( $arg{-v} ) ){ die "$version\n\n" }

## INPUT FILES/DATABASES
if( exists $arg{-d} ){ @filedb = split( /,/, $arg{-d} );}  ## input DBs
if( exists $arg{-i} ) { 
	@filein = split( /,/, $arg{-i} );
}
 
if( (scalar(@filein) + scalar(@filedb) ) < 1){ die "No input!!!\n",Help() }

## OUTPUT FILES

if( exists $arg{-o} and exists $arg{-text} ){
	open( OUT, ">$arg{-o}.txt" );
} elsif (exists $arg{-o}){
	open( OUT, ">$arg{-o}.htm" );
} else {
	open( OUT, ">&STDOUT" );
	
}

## My link styles

my $link_style =<<END;
<!--
a:link {text-decoration: none}
a:visited {text-decoration: none}
a:link {color: navy}
a:visited {color: teal}
-->
END

## start HTML
if(exists $arg{-text} ){
	print OUT "Domainer v4.0\n******************\n\n";
} else {
	print OUT start_html(-style=>{-code=>$link_style} );
	Header();
}

## Default variables

if( $arg{-E} =~ m/\d+/ ){ $e_val = $arg{-E} }
if( $arg{-T} =~ m/T/ ){ $total_score = "trusted_cutoff" }



################################# Main program #################################

## Creates a temporal file with DB's names

if( $arg{-d} ){
	if( !$arg{-u} or !$arg{-p} ){
		die "User name or password missing...!\n\n";
	}
    my $counter = 0;	
    foreach my $database ( @filedb ){
    	if( $arg{-I} =~ m/\d/ ){ $is_public = $arg{-I} } else {	$is_public = 1; }
	if( $database =~ m/(\S+):(\d)/ ){ ( $database, $is_public ) = ( $1, $2 ); }
	$db{ $database } = "above $total_score where is_public = $is_public";
	$filedb[$counter] = $database;
	$counter++;
	
	open( my $ofh, ">$database" );
	print STDERR "querying DB $database -- is_public = $is_public -- Total score above $total_score\n";
	
	if(exists $arg{-text} ){
		print OUT "Querying DB $database -- is_public = $is_public -- Total score above $total_score\n";
	} else {
		print OUT a("Querying DB $database -- is_public = $is_public -- Total score above $total_score"),br;
	}
	
	my $query ="
		SELECT e.accession, af.feat_name, eg.hmm_com_name
		FROM evidence e, egad..hmm2 eg, clone_info ci, asm_feature af, phys_ev pe
		WHERE af.feat_name = e.feat_name
		AND ci.asmbl_id = af.asmbl_id
		AND e.method = 'hmm2_search.dbi'
		AND ci.asmbl_id = af.asmbl_id
		AND eg.hmm_acc = e.accession
		AND pe.feat_name = af.feat_name
		AND pe.ev_type = 'working'
		AND ci.is_public = $is_public
		AND convert(float,e.expect_domain) < $e_val
		AND convert(float,e.total_score) > eg.$total_score";

	my $query_prok ="
		SELECT e.accession, e.feat_name, eg.hmm_com_name
		FROM feat_score f1, feat_score f2, evidence e, egad..hmm2 eg
		WHERE f1.score_id = 143
		AND f1.input_id = e.id
		AND f2.score_id = 144
		AND f2.input_id = e.id
		AND e.ev_type = 'HMM2'
		AND convert(float, f1.score) >= eg.$total_score
		AND eg.hmm_acc=e.accession
		AND e.feat_name like 'ORF%'";


	if(exists $arg{-prok}){
		$query = $query_prok;
	}	

		
	## SQL QUERY

	$dbh = DBI -> connect( "dbi:Sybase:server=SYBTIGR; packetSize=8092", $arg{-u}, $arg{-p});
	$dbh -> do( "use $database" ) || die "**** Cannot open db: $database ****\n";
	my $sth = $dbh -> prepare( $query );
	
        ## Do SQL QUERY

	$sth -> execute();
	$sth -> bind_columns( \$accession, \$feat_name, \$hmm_com_name );

	## Save temp DB file

	while( $sth -> fetch() ){
	    my @row=( $accession, $hmm_com_name, $feat_name );
	    print $ofh join( "\t", @row ),"\n";
	}
	
	# Close DB
	
	$sth -> finish;
	$dbh -> disconnect;
	close $ofh;
   }
}

## Load Pfam domains

my $cube; 		## HoHoH that stores pfam_ids, species_numbers and gene feat_names as keys and number of domains per pfma per gene per species as value
my %pfam_description; 	## link pfam_ids with pfam descriptions
my $sp_number = 0;	
my %pfam_groups = [];	## HoA that links pfam_ids as keys with species numbers as values in the array 
my %species_names = ();	## Link species numbers as keys with species names as values
@filein = (@filedb,@filein); ## input files

foreach my $species (@filein){ 
	open( IN, "<$species" )|| die "I cannot find file $species: $!\n";
	$sp_number++;
	$species_names{ $sp_number } = $species;
	
	while(<IN>){
		if(exists $arg{-notigr} ){next if m/TIGR\d{5}/;} ## Do not include TIGRFAMs
	     	chomp;
		next if(m/^[\s*\t*]$/);
	     	my ( $pfam, $com_name, $gene ) = split /\t+/;
		$pfam =~ s/\s+//g;
	     	$gene =~ s/\s+//g;
		$cube->{$pfam}{$gene}{$sp_number}++;
		$pfam_description{$pfam} = $com_name;
		$pfam_groups{$pfam}->[$sp_number] = $sp_number;
		
	}
	close IN;
}


## Store pfam IDs by groups

my %groups;	## HoA that links group_types (i.e. 1+2+3) as keys with pfam_ids as values of the array  
foreach my $pfam (keys %pfam_groups){
	$pfam_groups{$pfam} = join('+', grep( m/^\d+$/, @{ $pfam_groups{$pfam} } ) ); 
	push @{ $groups{ $pfam_groups{$pfam} } }, $pfam;
	delete ($groups{''});
}

#################
#### Results ####
#################



## Species description and total number of domains per species

my $number_of_domains = countDomains($cube);
if(exists $arg{-text} ){
	print OUT '='x100,"\n";
	print OUT "Group descriptions:\n";
} else {
	print OUT p,hr, h3('Group descriptions:'),p;
}
my $counter = 0;
foreach my $species (@filein){
	$counter++;
	my $current_number_of_domains = $number_of_domains->{$counter};
	if(exists $arg{-text} ){
		print OUT "Species number $counter: $species\tTotal number of domains:\t$current_number_of_domains\n";
	} else {
		print OUT a("Species number $counter: $species ($current_number_of_domains domains)"),br;
	}
}



## Number of domains shared by group and number of genes per group per species

if(exists $arg{-text} ){
	print OUT '='x100,"\n";
	print OUT "Venn diagram\n\n";
} else {
	print OUT p, hr, h3("Venn diagram"), p;
}

foreach my $this_group ( sort{ length($b) <=> length($a) } keys %groups){
	my $group_length = scalar( @{ $groups{ $this_group } } );
	my @species = split( /\+/, $this_group);
	my @genes_per_sp;
	my @number_of_genes;
	foreach my $species (@species){ 
		foreach my $pfam ( @{ $groups{ $this_group } } ){
			foreach my $gene ( keys %{ $cube->{$pfam} } ){
				push @{$genes_per_sp[$species]}, $gene if( exists( $cube->{$pfam}{$gene}{$species} ) );
			}
		}
		@{ $genes_per_sp[$species] } = (deleteDuplicates(\@{$genes_per_sp[$species]}));
		$number_of_genes[$species] = scalar( @{ $genes_per_sp[$species] } );
		
		## asign internal reference to list of genes per group per species
		
		my $flag_ref = $this_group.'_'.$species.'_';
		$number_of_genes[$species] = a({-href=>'#'.$flag_ref}, $number_of_genes[$species]) if ($number_of_genes[$species] > 0 and !exists $arg{-text});
	}
	@number_of_genes = grep(/\d+/, @number_of_genes);
	my $header = 'Domains specific for species';
	if($this_group=~ m/\+/){$header = 'Domains shared by species';}
	if(exists $arg{-text} ){
		print OUT "$header $this_group\t$group_length\tNumber of genes:\t",join("\t", @number_of_genes),"\n";
	} else {
		print OUT "$header ",
			a({-href=>'#'.$this_group}, $this_group),
			": $group_length (",
			join(' genes,', @number_of_genes),
			" genes)",
		  	br;
	}
}

## Lists of pfam_id, Description, number of domains, number of genes


foreach my $this_group ( sort{ length($b) <=> length($a) } keys %groups){
	my $header = 'Domains specific for species ';
	if($this_group=~ m/\+/){$header = 'Domains shared by species';}
	if(exists $arg{-text} ){
		print OUT '='x100,"\n";
		print OUT "$header shared by\t$this_group\n";
	} else {
		print OUT p, hr;
		print OUT a({-name=>$this_group}, h3("$header $this_group\n") );
	}
	
	my @species = split( /\+/, $this_group);
	my $domain_ref;
	
	foreach my $pfam ( @{ $groups{ $this_group } } ){
		my @genes_per_sp;
		my @domains_per_species;
			
		foreach my $species (@species){ 
			foreach my $gene ( keys %{ $cube->{$pfam} } ){
				if( exists( $cube->{$pfam}{$gene}{$species} ) ){
					$genes_per_sp[$species]++;
					$domains_per_species[$species] += $cube->{$pfam}{$gene}{$species};
				}
			}
			
		}
		@genes_per_sp = grep(/\d+/, @genes_per_sp);
		@domains_per_species = grep(/\d+/, @domains_per_species);
		
		## print OUT "$pfam: $pfam_description{$pfam}. Number of domains:".join(',',@domains_per_species)." Number of genes:".join(',',@genes_per_sp),
			br;
		$domain_ref = "http://www.sanger.ac.uk/cgi-bin/Pfam/getacc?$pfam" if($pfam =~ /PF/);
		$domain_ref = "http://www.tigr.org/tigr-scripts/CMR2/hmm_report.spl?acc=$pfam&user=access&password=access" if($pfam =~ /TIGR/);
		if(exists $arg{-text} ){
			print OUT "$pfam\t$pfam_description{$pfam}\tNumber of domains:\t",
			join("\t",@domains_per_species),
			"\tNumber of genes:\t", join("\t",@genes_per_sp ),"\n";
		} else {
			print OUT a( {-href=>$domain_ref, -target=>'_new'}, "$pfam"),
			": $pfam_description{$pfam}. (# domains: ",
			join(',',@domains_per_species),")  (# genes: ",
			a({-href=>'#'.$pfam}, join(',',@genes_per_sp) ),")",
			br;
		}
	}
	
}

## List pfam_id, species and gene feat_names

 foreach my $this_group ( sort{ length($b) <=> length($a) } keys %groups){
 	print OUT p,h3("Genes with domains shared by $this_group") unless (exists $arg{-text} );
	
	my @species = split( /\+/, $this_group);
	
	foreach my $pfam ( @{ $groups{ $this_group } } ){
		if(exists $arg{-text} ){
			# print OUT '='x100,"\n";
		} else {
			print OUT hr; 
		}
		foreach my $species (@species){
			my @genes_per_sp;
			my $domain_ref;
			foreach my $gene ( keys %{ $cube->{$pfam} } ){
				if( exists( $cube->{$pfam}{$gene}{$species} ) ){
					my $feat_name_sybil = getSybilfeat_name($gene, $species_names{$species});
					my $prok_manatee = a({-target=>'_new', -href=>"http://manatee.tigr.org/tigr-scripts/prok_manatee/shared/ORF_infopage.cgi?db=$species_names{$species}&orf=$gene"},$gene);
					my $manatee_link = a({-target=>'_new', -href=>"http://manatee.tigr.org/tigr-scripts/euk_manatee/shared/ORF_infopage.cgi?db=$species_names{$species}&orf=$gene"},$gene);
					my $sybil_link = a({-target=>'_new',-class=>'button', -href=>"http://manatee.tigr.org/tigr-scripts/sybil-apx2/showProtein.pl?db=apx2&uniquename=$feat_name_sybil"},"<small> s </small>"); 
					$manatee_link = $prok_manatee if ( $arg{-prok} eq 'T');
					push @genes_per_sp, $manatee_link . $sybil_link;
				}
			}
			$domain_ref = "http://www.sanger.ac.uk/cgi-bin/Pfam/getacc?$pfam" if($pfam =~ /PF/);
			$domain_ref = "http://www.tigr.org/tigr-scripts/CMR2/hmm_report.spl?acc=$pfam&user=access&password=access" if($pfam =~ /TIGR/);
				
			print OUT a({-name=>$pfam, -target=>'_new', -href=> $domain_ref}, $pfam),
				" ($species_names{$species}): ",
				join(', ',@genes_per_sp),
				p unless (exists $arg{-text} );
			
		}
	}
}

## List all gene feat_names grouped by group and species

foreach my $this_group ( sort{ length($b) <=> length($a) } keys %groups){
	my @species = split( /\+/, $this_group);
	my @genes_per_group_per_sp;
	my @number_of_genes;
	print OUT hr unless (exists $arg{-text} );
	foreach my $species (@species){ 
		foreach my $pfam ( @{ $groups{ $this_group } } ){
			foreach my $gene ( keys %{ $cube->{$pfam} } ){
				push @{$genes_per_group_per_sp[$species]}, $gene if( exists( $cube->{$pfam}{$gene}{$species} ) );
			}
		}
		@{ $genes_per_group_per_sp[$species] } = (deleteDuplicates(\@{$genes_per_group_per_sp[$species]}));
		my $flag_ref = $this_group.'_'.$species.'_';
		print OUT a({-name=>$flag_ref}, h3("Genes from species $species ($species_names{$species}) with domains shared by species $this_group:")),
		p unless (exists $arg{-text} );
		print OUT join(', ',@{$genes_per_group_per_sp[$species]}),
		p unless (exists $arg{-text} );
		print OUT hr unless (exists $arg{-text} );

	}
	
	
}
print OUT end_html unless (exists $arg{-text} );
close OUT;

## clean temp files

unless(exists $arg{-keepdb}){
	my $rest;
	foreach my $file (@filedb){
		my $cmd = "rm $file";
		my $err = system($cmd);
		print STDERR "$err\n\n" if $err;
	}	
}

################
## Subrutines ##
################

## Get a reference to an array with duplicated elements and return single elements by order of aperance

sub deleteDuplicates {
	my ($array_of_duplicated_elements) = @_;
	my %hash_of_single_elements;
	my $counter = 0;
	foreach ( @{ $array_of_duplicated_elements} ){
		s/\s+//g;
		$hash_of_single_elements{$_} = $counter;
		$counter++;
	}
	sort{ $hash_of_single_elements{$a} <=> $hash_of_single_elements{$b} } keys %hash_of_single_elements
}

sub getSybilfeat_name {
	my ($feat_name, $db) = @_;
	my ($asmbl, $gene) = split ('.m', $feat_name);
	return $db.'_'.$asmbl.'_'.$asmbl.'.p'.$gene;
}

sub countDomains {
	my ($cube) = @_;
	my $duplicate;
	my %number_of_domains;
	foreach my $pfam (keys %{$cube}){
		foreach my $gene (keys %{ $cube->{$pfam} } ){
			foreach my $sp_number (keys %{$cube->{$pfam}{$gene} } ){
				$number_of_domains{$sp_number}++ unless exists $duplicate->{$pfam}{$sp_number};			
				$duplicate->{$pfam}{$sp_number} = 1;
			}
		}
	}
	return \%number_of_domains;
	
}

## Print header

sub Header{
my $header =<<END;
Domainer v4.0<br>
****************<br>	
<br>
END
 
print OUT a($header), 
}

## Print help

sub Help{
print <<START


Description: 
-----------

	Domainer determines the number of domains shared by two or more species.
	

	Usage: 
			-From flat files:
			
			$0 -i <inputfile1,inputfile2,etc> [options]

			-From DBs:
			
			$0 -d <db1,db2,etc> [options]

	Input files: list of pfam_ids, gene_ids, pfams_descriptions (in columns, tab delimited) for each species.
	
	Options:
	-------
 	
        -text T : generate a tab delimited text file
	
	-o output file [STDOUT]
	
		
	DATABASE Options:
	
        -d query pfam domains at a local Sybase database
	
		-d pva1,pfa1,pya1    -> Search pva1, pfa1 and pya databases
		It is also possible to set is_public values for each database adding that value after the DB name separated by ':':
		-d pva1,pfa1:0,pya1  -> Search pfa1 database using is_public = 0 and the default values for pva1 and pya1
	
	-p password
	
	-u username
	
	-I set is_public value for all DB searched [1]
	
	-E set e-value of hmm hit [0.001]
	
	-T [T/F] show domains above trusted cutoff [F:gathering cutoff]

	-notigr [T/F] do not include TIGRFAMs entries in the analysis [F]	
	
	-keepdb [T/F] do not remove temp file with database query output [F]
	
	-v print version
	
	-prok [T/F] query a prokaryote database [F]	

START
}

END {$dbh -> disconnect if($dbh)}
