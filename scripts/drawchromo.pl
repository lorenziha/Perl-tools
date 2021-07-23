#!/usr/local/bin/perl
use GD;

$help = "Usage: $0 -i coordinates_file -d [all , list domain1:color1:shape1,domain2:color2:shape2..domainN:colorN:shapeN] [-W # -H # -s # -scale # -n # -c dom1,dom2,dom3,dom4,dom5 -g picture_name] > file_name.png\n
	\ndrawchromo.pl draws a series of chromosomes using as an imput a txt file containing the coordinates for each domain:
	The format of the txt file is as follows:\n
	>Chromosome_ID, length , [comments]
	domain1_start , domain1_end , Domain ID , color , shape, motif/domain seq/comments (they will appear on the references)\n
	Example:\n
	>Chromosome:1,1595673
	7234,10026,gene1,a,r,kinase
	428763,455507,gene2,a,r,HSP70
	>Chromosome:2,3098763
	227222,249893,gene3,a,r,GAPDH\n
	Color values [a,1-350]:
	a -> drawprot.pl asign the domain color automaticaly\n
	Shapes values [r,c, and more...]:
	r -> draw domains as rectangles
	c -> draw domains as ellipses\n
	Optional parameters:\n	
		-W # picture width [1500]
		-H # picture high [1500]
		-d list of domains to be draw. all -> all domains
		   It is possible to set colors and shapes for every domain using ':' as separators
		   So far, the shapes supported are: r,c,p1,p2,p3,p4 and p5
		-c It centers the proteins based on domain 1,2,3,4 or 5 whatever it comes first.
		-g open an old picture (png) file.
		-scale [0..1] scales the picture. 1 is the original size.
		-s # set the distance between proteins [15]
		-n # shows all chromosomes with a particular length\n";
			


die $help if(grep(/-h/,@ARGV));

## create a new image

my %arg = @ARGV;
my $page_width = 1500;
my $page_high = 1500;
my %center;
my $motif_num = 0;
my $margin = 180;
my $previous = 'long';
my $max = 0;
 
if($arg{-W} =~ /\d+/){
	$page_width = $arg{-W};
}
if($arg{-H} =~ /\d+/){
	$page_high = $arg{-H};
}
$im = new GD::Image($page_width,$page_high);


if($arg{-d} =~ /\S+/){
	my @domains = split(/,/,$arg{-d});
	foreach my $key (@domains){
		($name,$color,$shapes) = split (/:/,$key);
		$domains{$name} = $color;
		$domains{$name} = "C" if ($color !~ /\d+/);
		$shape{$name} = $shapes if $shapes =~ m/\S+/;
		
	} 
}
my $step = $arg{-s} || 15;
my $scale = $arg{-scale} || 1;

## allocate some colors

open (COLOR,"</usr/X11/lib/X11/rgb.txt") || die "I cannot locate RGB table\n";
while(<COLOR>){
	chomp;
	$color_num++;
	s/^\s+//g;
	($red,$green,$blue,$color_id) = split /\s+|\t+/;
	if ($color_num == 2){
		$color_num = 0;
		$color{$color_num} = $im->colorAllocate($red,$green,$blue) unless ($blue !~ /\d+/);
		push @color, $color{$color_num};
	}
}	
$black = $color[20];
$while = $color[1];

close (COLOR);

# make the background interlaced

$im->interlaced('true');

## Load dwp file

open(my $fhi,"<".$arg{-i}) || die "$help\n";
while(<$fhi>){
	chomp;
	my $line++;
	
	## Chromosome ID
		
	if(/^>/){
		$counter = Max($counter,$max);
		$counter += $step;
		@gene = split/,/;
		$header = shift @gene;
		$header =~ s/>//;
		my $chromosome_len = shift @gene;
		$num_motifs = shift @gene;
		$eval = shift @gene;
		
		## Normalize chromosome length (-n option)
		
		if (exists $arg{-n}){
			my $normalization_val = $arg{-n};
			$scale = $normalization_val/$chromosome_len;
		}
		$chromosome_len *= $scale;
		
		## Draw center line
		
		$im->filledRectangle($margin,$counter,$chromosome_len+$margin, $counter,$black);
		
		## Add telomeres
		
		$im->setThickness(2);
		telomere_r(\$im,$margin,$chromosome_len+$margin,$counter,$black) if exists $arg{-t};
		telomere_l(\$im,$margin,$chromosome_len+$margin,$counter,$black) if exists $arg{-t};
		$im->setThickness(1);
				
		## Write Labels
		
		$im->string(gdGiantFont,10,$counter-5,$header,$black) unless ($arg{-g} =~ /\S+/);
		
	}
	
	## Gene coordinates
	
	elsif(/^\d+,\d+,\S+,\S+,\w/){
		chomp;
		@motif = split/,/;
		$start = $scale * (shift @motif) + $margin;
		$end   = $scale * (shift @motif) + $margin;
		$motif_id = shift @motif;
		$color_id = shift @motif;
		if( exists($shape{$motif_id}) ){
			$shape = $shape{$motif_id};
			shift @motif;
		}
		else {
			$shape = shift @motif;
		}
		$seq = shift @motif;
		$motif_seq{$motif_id} = $seq;
		$strand = shift @motif || 1;
		
		## Color asignment: 
				
		if($color_id =~ /a/){
			if(!$motif_color{$motif_id}){
				$motif_num++;
				if($domains{$motif_id} =~ /\d/){
					$motif_color{$motif_id} = $domains{$motif_id};
				}
				else {
					$motif_color{$motif_id} = $color[$motif_num*10];
				}
				
			}
		}
		else {
			$motif_color{$motif_id} = $color[$color_id];
		}

		
		$motif_shape{$motif_id} = $shape;
		if($shape =~ /r$/ and ($domains{$motif_id} or $arg{-d} =~ /all/)){
			
			$im->filledRectangle($start , $counter-5 , $end , $counter+5 , $motif_color{$motif_id});
			$im->rectangle($start , $counter-5 , $end , $counter+5 , $black);

		}
		elsif($shape =~ /r2/){
			$im->setThickness(2);
			## Move legend downward if block is too small
			
			my $count_bkp = $counter;
			my $midpoint = $start+int(($end-$start)/2);
			for (my $i = -55; $i < 55 ; $i += 1){
				my $pixel =  $im->getPixel($i+$midpoint,$count_bkp+14);
				my ($r, $g, $b) = $im->rgb($pixel);
				if ($r == 0 and $g == 0 and $b == 0){
					while($r == 0 and $g == 0 and $b == 0){
						$count_bkp += 6;
						$pixel =  $im->getPixel($i+$midpoint,$count_bkp+14);
						($r, $g, $b) = $im->rgb($pixel);
					}
					$i = -55;
					
				}
				$max = Max($max,$count_bkp);
			}
			
			
			r2(\$im,$start,$end,$counter,$count_bkp,$motif_color{$motif_id},$seq,$strand);
			$im->setThickness(1);
		}
		elsif($shape =~ /rna/){
			if( Module($end - $start) < 5){
				$end = $start + 5 if $start < $end;
				$end = $start - 5 if $start >= $end;
			}
			rna(\$im,$start,$end,$counter,$motif_color{$motif_id},$seq);
		}
		elsif($shape =~ /c/){
			$width = $end - $start;
			$centerx = ($width/2)+$start;
			$im->filledEllipse($centerx,$counter,$width,10,$motif_color{$motif_id}) if($domains{$motif_id} or $arg{-d} =~ /all/); 
			$im->ellipse($centerx,$counter,$width,10,$black) if($domains{$motif_id} or $arg{-d} =~ /all/);
			
		}
		elsif($shape =~ /p1/){
			pin1(\$im,$start,$end,$counter,$motif_color{$motif_id});
		
		}		
		elsif($shape =~ /p2/){
			pin2(\$im,$start,$end,$counter,$motif_color{$motif_id});
		
		}		
		elsif($shape =~ /p3/){
			$im->setThickness(2);
			pin3(\$im,$start,$end,$counter,$motif_color{$motif_id});
			$im->setThickness(1);
		}		
		elsif($shape =~ /p4/){
			pin4(\$im,$start,$end,$counter,$motif_color{$motif_id});
		
		}		
		elsif($shape =~ /p5/){
			pin5(\$im,$start,$end,$counter,$motif_color{$motif_id});
		
		}		
		elsif($shape =~ /line/){
			line(\$im,$start,$end,$counter,$motif_color{$motif_id});
		
		}		
		
		
	}
	
		
}

## Draw references

$counter += $step;

foreach $key (sort {$a <=> $b} keys %motif_color){
	$counter += 25 if($domains{$key} or $arg{-d} =~ /all/);
	if($motif_shape{$key} =~/^r$/){ 
		$im->filledRectangle($margin+$adjust , $counter-5 ,191+$adjust , $counter+5 , $motif_color{$key} ) if($domains{$key} or $arg{-d} =~ /all/); 
		$im->rectangle($margin+$adjust , $counter-5 ,191+$adjust , $counter+5 , $black ) if($domains{$key} or $arg{-d} =~ /all/);
	}
	elsif($motif_shape{$key} =~/c/) { 
		$im->filledEllipse($margin+$adjust+5,$counter,15,10,$motif_color{$key}) if($domains{$key} or $arg{-d} =~ /all/); 
		$im->ellipse($margin+$adjust+5,$counter,15,10,$black) if($domains{$key} or $arg{-d} =~ /all/);
	}
	elsif($motif_shape{$key} =~/p1/){pin1(\$im, $margin+$adjust, 1+$margin, $counter+10,$motif_color{$key})}
	elsif($motif_shape{$key} =~/p2/){pin2(\$im, $margin+$adjust, 1+$margin, $counter-5,$motif_color{$key})}
	elsif($motif_shape{$key} =~/p3/){pin3(\$im, $margin+$adjust, 1+$margin, $counter,$motif_color{$key})}
	elsif($motif_shape{$key} =~/p4/){pin4(\$im, $margin+$adjust, 1+$margin, $counter+5,$motif_color{$key})}
	elsif($motif_shape{$key} =~/p5/){pin5(\$im, $margin+$adjust, 1+$margin, $counter-5,$motif_color{$key})}
	elsif($motif_shape{$key} =~/line/){line(\$im, $margin+$adjust, 1+$margin, $counter,$motif_color{$key})}
	elsif($motif_shape{$key} =~/rna/){rna(\$im, $margin+$adjust-5, $margin+10+$adjust, $counter+5,$motif_color{$key})}
		
	$im->string(gdGiantFont,198+$adjust,$counter-5,"$key- Family ".$motif_seq{$key},$black) if($domains{$key} or $arg{-d} =~ /all/);
}

## make sure we are writing to a binary stream
binmode STDOUT;

## Convert the image to PNG and print it on standard output
print $im->png;	

###########
## Shapes #
###########

sub Max{
	my ($max, $count_bkp) = @_;
	if($max >= $count_bkp){return $max;}
	else {return $count_bkp;}
}
sub pin1{
	my ($im,$start,$end,$counter,$color) = @_;
	my $width = $end - $start;
	my $centerx = ($width/2)+$start;
	${$im}->filledEllipse($centerx,$counter-15,7,7,$color); 
	${$im}->line($centerx,$counter,$centerx,$counter-10,$color);

}
sub pin2{
	my ($im,$start,$end,$counter,$color) = @_;
	my $width = $end - $start;
	my $centerx = ($width/2)+$start;
	${$im}->filledEllipse($centerx,$counter+15,7,7,$color); 
	${$im}->line($centerx,$counter,$centerx,$counter+10,$color);

}
sub pin3{
	my ($im,$start,$end,$counter,$color) = @_;
	my $width = $end - $start;
	my $centerx = ($width/2)+$start;
	${$im}->filledEllipse($centerx,$counter,13,13,$white); 
	${$im}->ellipse($centerx,$counter,13,13,$color);
}
sub pin4{
	my ($im,$start,$end,$counter,$color) = @_;
	my $width = $end - $start;
	my $centerx = ($width/2)+$start;
	${$im}->filledRectangle($centerx-3,$counter-11,$centerx+3,$counter-5,$color);
	${$im}->line($centerx,$counter,$centerx,$counter-5,$color);
 
}
sub pin5{
	my ($im,$start,$end,$counter,$color) = @_;
	my $width = $end - $start;
	my $centerx = ($width/2)+$start;
	${$im}->filledRectangle($centerx-3,$counter+5,$centerx+3,$counter+11,$color);
	${$im}->line($centerx,$counter,$centerx,$counter+5,$color);
}
sub line{
	my ($im,$start,$end,$counter,$color) = @_;
	my $width = $end - $start;
	my $centerx = ($width/2)+$start;
	${$im}->line($centerx,$counter-3,$centerx,$counter+3,$color);
}

## Labeled rectangle with arrows

sub r2{
	my ($im,$start,$end,$counter,$count_bkp,$color,$seq,$strand) = @_;
	my $mid = int(($start+$end)/2);
	my $len = 4 * length($seq); 
	${$im}->filledRectangle($start, $counter-8 , $end, $counter+8 , $color);
	${$im}->rectangle($start, $counter-8 , $end , $counter+8 , $black);
	${$im}->string(gdLargeFont,$mid-$len,$count_bkp+8,$seq,$black);
	
	## arrows
	
	if($strand < 0){
		${$im}->line($mid-10,$counter-13,$mid+10,$counter-13,$black);
		${$im}->line($mid-10,$counter-13,$mid-7,$counter-16,$black);
		${$im}->line($mid-10,$counter-13,$mid-7,$counter-10,$black);
	} elsif($strand > 0){
		${$im}->line($mid-10,$counter-13,$mid+10,$counter-13,$black);
		${$im}->line($mid+7,$counter-16,$mid+10,$counter-13,$black);
		${$im}->line($mid+7,$counter-10,$mid+10,$counter-13,$black);
	}
}

## RNA genes

sub rna{
	my ($im,$start,$end,$counter,$color,$seq) = @_;
	my $strand = 1;
	if($end < $start){$strand = -1}
	if($strand > 0){
		${$im}->filledRectangle($start, $counter - 10 , $end, $counter, $color);
		${$im}->line($end-5,$counter-14,$end,$counter-10,$color);
	} else {
		${$im}->filledRectangle($end, $counter , $start, $counter+10, $color);
		${$im}->line($end+5,$counter+14,$end,$counter+10,$color);
	}
}

sub telomere_r {
	my ($im,$start,$end,$counter,$color) = @_;
	my $ul_point = $counter - 8;
	my $ll_point = $counter + 8;
	my $mr_point = $counter;
	my $r_telomere = new GD::Polygon;
	${$im}->line($end,$ul_point,$end,$ll_point,$color);
	${$im}->line($end,$ul_point,$end+20,$mr_point,$color);
	${$im}->line($end+20,$mr_point,$end,$ll_point,$color);
}
sub telomere_l {
	my ($im,$start,$end,$counter,$color) = @_;
	my $ul_point = $counter - 8;
	my $ll_point = $counter + 8;
	my $mr_point = $counter;
	my $r_telomere = new GD::Polygon;
	${$im}->line($start,$ul_point,$start,$ll_point,$color);
	${$im}->line($start,$ul_point,$start-20,$mr_point,$color);
	${$im}->line($start-20,$mr_point,$start,$ll_point,$color);
	
}
sub Module{
	$_[0] *= -1 if $_[0] < 0;
	return $_[0];
}
