#!/usr/local/bin/perl 
use lib ($ENV{EUK_MODULES}, $ENV{EGC_SCRIPTS});
use Egc_library;
use strict;
use DBI;
use Gene_obj;
use Fasta_reader;
use Data::Dumper;

my $usage = "usage: $0 -d db -a asmbl_accessions -verbose T/F [F] -p T/F [F, just_print_protein_seqs]\n\n";
my %arg = @ARGV;
my $db = uc($arg{-d}) or die $usage;
my $asmbl_accessions = $arg{-a} or die $usage;
my $VERBOSE = $arg{-verbose} || 'F';
my $PROTEIN = $arg{-p} || 'F';



#############################################################################
# begin program

my $current_dir = `pwd`;
chomp $current_dir; 

open(ASMBL,"<$asmbl_accessions") || die "ERROR, I cannot open $asmbl_accessions: $!\n\n";
while(<ASMBL>){
	chomp;
	next if m/^[\s\t]*$/;
	my ($asmbl_id, $query, $strand) = split m/\t+/;
	$asmbl_id =~ s/\s+//g;
	my ($legacy_asmbl_id) = $asmbl_id =~ m/_(\d+)_/; 
	my $sequence_dir = "/usr/local/annotation/$db/asmbls/$legacy_asmbl_id";
	my $genewise_dir = $current_dir."/TARGET_GENOME/$asmbl_id";
	my $genewise_file = `ls $genewise_dir/$asmbl_id*.gff3`;
	chomp $genewise_file;
	my $sequence = $sequence_dir."/$legacy_asmbl_id.contig";
	unless (-d 	$sequence_dir){
		warn "ERROR, I cannot find $sequence_dir\n";
		next;
	}
	unless (-d 	$genewise_dir){
		warn "ERROR, I cannot find $genewise_dir\n";
		next;
	}
	unless (-s 	$genewise_file){
		warn "ERROR, I cannot find $genewise_file\n";
		next;
	}
	unless (-s 	$sequence){
		warn "ERROR, I cannot find $sequence\n";
		next;
	}
	
	## Get asmbl start coordinate
	my ($asmbl_part) = $asmbl_id =~ m/_(\d+)$/;
	#print "ASMBL_PART= $asmbl_part\n";
	unless ($asmbl_part){
		$asmbl_part = 1;
	}
	#print "ASMBL_PART= $asmbl_part\n";
	my $asmbl_start_coord = ($asmbl_part - 1) * 880;
	#print "ASMBL_PART= $asmbl_start_coord\n";
	## Get asmbly sequence
	my $fasta_reader = new Fasta_reader($sequence);
	my $seq_obj = $fasta_reader->next();
	my $asmbl_sequence = $seq_obj->get_sequence();
	
	#print "REF= ",ref( $asmbl_sequence),"\n",length($asmbl_sequence),"\n";
	
	## create gene_obj from predictions
	open(GENEWISE,"<$genewise_file") || die "I cannot open $genewise_file: $!\n\n";
	my @genewise = (<GENEWISE>);
	close GENEWISE;
	
	#print "genewise= ",join('-',@genewise),"\n";
	
	my @gene_objs = &get_gene_obj_from_genewise(\$asmbl_sequence,$asmbl_start_coord, @genewise);
	foreach my $gene_obj (@gene_objs){
		$gene_obj->{asmbl_id} = $legacy_asmbl_id;
		if($VERBOSE eq 'T'){
			print 'TU_feat_name= ',$gene_obj->{TU_feat_name},"\n";
			print 'Model_feat_name= ',$gene_obj->{Model_feat_name},"\n";
			print 'protein_seq= ',$gene_obj->{protein_seq},"\n";
			print 'CDS_sequence= ',$gene_obj->{CDS_sequence},"\n";
			print 'gene_span= ',"@{$gene_obj->{gene_span}}","\n";
			print 'model_span= ',"@{$gene_obj->{model_span}}","\n";
			print 'num_exons= ',$gene_obj->{num_exons},"\n";
			print 'is_5prime_partial= ',$gene_obj->{is_5prime_partial},"\n";
			print 'is_3prime_partial= ',$gene_obj->{is_3prime_partial},"\n";
			print 'strand= ',$gene_obj->{strand},"\n\n";
		}
		unless ($PROTEIN eq 'T'){
			
			eval { $gene_obj->set_CDS_phases(\$asmbl_sequence) };
			if ($@){
				warn "\nWARNING, error getting phases from $gene_obj->{TU_feat_name}: $@\n";
			}
			else {
				print $gene_obj->to_GFF3_format();
				print "\n";
			}
		}
		else {
			print ">$gene_obj->{TU_feat_name} $gene_obj->{asmbl_id} @{$gene_obj->{gene_span}} 5_part $gene_obj->{is_5prime_partial} 3_part $gene_obj->{is_3prime_partial} \n";
			print "$gene_obj->{protein_seq}\n";
		}
		
	}
}
close ASMBL;
exit(0);

#############################################################
sub get_gene_obj_from_genewise {
	my ($sequence_ref, $asmbl_start_coord, @genewise) = @_;
	my (@model, @model_obj);
	my ($start_flag, $stop_flag, $end5_partial, $end3_partial, $end5, $end3, $exon_count) = (0,0,0,0,0,0,0);
	my $line_counter = 0;
	my ($old_gene_id, $old_transcript_id);
	my $source = '';
	foreach my $line (@genewise){
		chomp $line;
		$line_counter++;
		
		next if ($line =~ m/^#/);
		next if ($line =~ m/^[\s\t]*$/);
		#print "LINE= $line\n";
		my $feat = &get_genewise_feature($line);
		$source = $feat->{analysis} || $source;
		#print "DUMPER= ",Dumper($feat);
		
		if ($old_gene_id && ($feat->{gene_id} ne $old_gene_id || $line =~ m/^\//)){
			## NEW GENE
			
			my $gene_obj = &create_gene_obj($sequence_ref, $asmbl_start_coord, @model);
			($end5_partial, $end3_partial) = &get_5_3_partial($gene_obj->{protein_seq});
			
			$gene_obj->{TU_feat_name} = $old_gene_id;
			$gene_obj->{Model_feat_name} = $old_transcript_id;
			$gene_obj->{is_5prime_partial} = $end5_partial;
			$gene_obj->{is_3prime_partial} = $end3_partial;
			$gene_obj->{source} = $source;
			
			#print Dumper($gene_obj);
			push @model_obj, $gene_obj;
			## reset @model
			@model = ();
			
			## reset flags
			($start_flag, $stop_flag, $end5_partial, $end3_partial, $end5, $end3, $exon_count) = (0,0,0,0,0,0,0);
			
		}
		
		
		$old_gene_id = $feat->{gene_id};
		$old_transcript_id = $feat->{transcript_id};

		
		if ( $feat->{strand} eq '+' ){
			## Forward strand
			if ($feat->{feature} eq 'cds'){
				$exon_count++;
				$end5 = $feat->{lend};
				$end3 = $feat->{rend};
				push @model, { end5 => $end5 , end3 => $end3};
			}
		}
		else {
			## Reverse strand
			if ($feat->{feature} eq 'cds'){
				$exon_count++;
				$end5 = $feat->{lend};
				$end3 = $feat->{rend};
				push @model, { end5 => $end5 , end3 => $end3};
			}
		}
	}
	return @model_obj;
}

sub get_5_3_partial {
	my $protein_seq = @_;
	my $end5_partial = $protein_seq =~ m/^M/? 0 : 1;
	my $end3_partial = $protein_seq =~ m/\*$/? 0 : 1;
	return ($end5_partial, $end3_partial);
}

sub create_gene_obj {
	my ($sequence_ref, $asmbl_start_coord, @model) = @_;
	my (%mrna, %CDS);
	
	#print "asmbl_start_coord $asmbl_start_coord\n";
	foreach my $coord_pair (@model){
		#print Dumper($coord_pair);
		my ($end5, $end3) = ( $coord_pair->{end5}+$asmbl_start_coord,$coord_pair->{end3}+$asmbl_start_coord );
		$mrna{ $end5 } = $end3;
		$CDS { $end5 } = $end3;
		#print 'COORDS='.join('-',@{$coord_pair}),"\n";
	}
	my $gene_obj = new Gene_obj();
	$gene_obj->{com_name} = 'genewise prediction';
	$gene_obj->populate_gene_obj(\%CDS, \%mrna, $sequence_ref);
	return $gene_obj;
}

sub get_genewise_feature {
	my $line = shift;
	my @line = split(m/\t/,$line);
	
	return {	analysis => $line[1],
				feature => $line[2],
				asmbl	=> $line[0],
				lend	=> $line[3],
				rend	=> $line[4],
				score	=> $line[5],
				strand	=> $line[6],
				phase	=> $line[7],
				gene_id => $line[8],
				transcript_id => $line[8]
			};
}

