#!/usr/local/bin/perl
use GD;

$help = "drawprot.pl draw a series of proteins using as an imput a txt file containing the coordinates for each domain:
	The format of the txt file is as follows:\n
	>Protein ID , length , [comments]
	domain1_start , domain1_end , Domain ID , color , shape, motif/domain seq/comments (they will appear on the references)
	Example:\n
	>Pv001100|35.m00825|52371,159,6,1.38e-53
	7,26,5,a,r,LPFGCTEWDLSKQLTEWQI
	42,57,3,a,r,YIIWYYVHNYLRRKY
	>Pv083565|393.m00243|9913,281,1,5.70e-04
	227,249,1,a,r,MNRLWKWCHDLATKCGIPHEYK\n
	Color values [a,1-350]:
	a -> drawprot.pl asign the domain color automaticaly\n
	Shapes values [r,c]:
	r -> draw domains as rectangles
	c -> draw domains as ellipses\n	
	usage: drawplot -i coordinates_file -d [all , list domain1:color1:shape1,domain2:color2:shape2..domainN:colorN:shapeN] [-W # -H # -s # -c dom1,dom2,dom3,dom4,dom5 -g picture_name] > file_name.png\n
		-W # picture width [1500]
		-H # picture high [1500]
		-d list of domains to be draw. all -> all domains
		   It is possible to set colors and shapes for every domain using ':' as separators
		   So far, the shapes supported are: r,c,p1,p2,p3,p4 and p5
		-c It centers the proteins based on domain 1,2,3,4 or 5 whatever it comes first.
		-g open an old picture (png) file.
		-s # set the distance between proteins [15]\n";
	


die $help if(grep(/-h/,@ARGV));

## create a new image

my %arg = @ARGV;
my $page_width = 1500;
my $page_high = 1500;
my %center;

if($arg{-W} =~ /\d+/){
	$page_width = $arg{-W};
}

if($arg{-H} =~ /\d+/){
	$page_high = $arg{-H};
}

if($arg{-g} =~ /\S+/){  ## Open existent image
	
	$im = new GD::Image($arg{-g});
	$motif_num = $arg{-s}; ## number of previous domains -> avoid repeat the same colors
	$adjust = 250;
}
else{
	
	$im = new GD::Image($page_width,$page_high);
	$motif_num = 0;
}

if($arg{-d} =~ /\S+/){
	@domains = split(/,/,$arg{-d});
	
	foreach $key (@domains){
		($name,$color,$shapes) = split (/:/,$key);
		$domains{$name} = $color;
		$domains{$name} = "C" if ($color !~ /\d+/);
		$shape{$name} = $shapes if $shapes =~ m/\S+/;
	} 
}
my $step = $arg{-s} || 15;

my $margin = 180;


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

close (COLOR);

$black = $color[20];

# make the background interlaced

$im->interlaced('true');



open(my $fhi,"<".$arg{-i});

if($arg{-c} =~ /\S+/){
	($dom1,$dom2,$dom3,$dom4,$dom5) = split (/,/,$arg{-c});
	
	while(<$fhi>){
		chomp;
		
		
		if(/^>/){ ## Header
			$model++;
			@gene = split/,/;
			$header = shift @gene;
			$header =~ s/>//;
			
			
		}
		elsif(/^\d+,\d+,$dom1,\S+,\w/ ){
			
			chomp;
			
			@motif = split/,/;
			
			$start = shift @motif;
			$center{$header} = int($page_width/2)-$start;
			$start1 = $start if ($model == 1);
		}
		elsif(/^\d+,\d+,$dom2,\S+,\w/ and $model == 1){
			chomp;
			
			@motif = split/,/;
			
			$start2 = shift @motif;
			
			
		}
		elsif(/^\d+,\d+,$dom3,\S+,\w/ and $model == 1){
			chomp;
			
			@motif = split/,/;
			
			$start3 = shift @motif;
			
			
		}
		elsif(/^\d+,\d+,$dom4,\S+,\w/ and $model == 1){
			chomp;
			
			@motif = split/,/;
			
			$start4 = shift @motif;
			
			
		}
		elsif(/^\d+,\d+,$dom5,\S+,\w/ and $model == 1){
			chomp;
			
			@motif = split/,/;
			
			$start5 = shift @motif;
			
			
		}
		elsif(/^\d+,\d+,$dom2,\S+,\w/){
			chomp;
			
			@motif = split/,/;
			
			$start = shift @motif;
			$center{$header} = int($page_width/2)-$start-($start1-$start2) unless $center{$header};
			
		}
		elsif(/^\d+,\d+,$dom3,\S+,\w/){
			chomp;
			
			@motif = split/,/;
			
			$start = shift @motif;
			$center{$header} = int($page_width/2)-$start-($start1-$start3) unless $center{$header};
			
		}
		elsif(/^\d+,\d+,$dom4,\S+,\w/){
			chomp;
			
			@motif = split/,/;
			
			$start = shift @motif;
			$center{$header} = int($page_width/2)-$start-($start1-$start4) unless $center{$header};
			
		}
		elsif(/^\d+,\d+,$dom5,\S+,\w/){
			chomp;
			
			@motif = split/,/;
			
			$start = shift @motif;
			$center{$header} = int($page_width/2)-$start-($start1-$start5) unless $center{$header};
			
		}
	}

}
close ($fhi);

open(my $fhi,"<".$arg{-i});

while(<$fhi>){
	chomp;
	my $line++;
	
	if(/^>/){ ## Header
		
		$counter += $step;
		@gene = split/,/;
		$header = shift @gene;
		$header =~ s/>//;
		if($arg{-c} =~ /\S+/){
			
			$center = $center{$header};
			
		}
		else {
			$center = 0;			
		}
		
		$gene_length = shift @gene;
		$num_motifs = shift @gene;
		$eval = shift @gene;
	
		$im->filledRectangle($margin+$center,$counter-1,$gene_length+$margin+$center, $counter+1,$black) unless ($arg{-g} =~ /\S+/);
		$im->string(gdSmallFont,10,$counter-5,$header,$black) unless ($arg{-g} =~ /\S+/);
		
	}
	elsif(/^\d+,\d+,\S+,\S+,\w/){
		chomp;
		
		@motif = split/,/;
			
		$start = shift @motif;
		$end = shift @motif;
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
				
		if($shape =~ /r/ and ($domains{$motif_id} or $arg{-d} =~ /all/)){
			#print STDERR "$motif_id:$motif_color{$motif_id}\n";
			$im->filledRectangle($margin+$start+$center , $counter-5 , $margin+$end+$center , $counter+5 , $motif_color{$motif_id});
			#$im->rectangle($margin+$start+$center , $counter-5 , $margin+$end+$center , $counter+5 , $black);

		}
		
		elsif($shape =~ /c/){
			
			$width = $end - $start;
			$centerx = ($width/2)+$start;
			$im->filledEllipse($margin+$centerx+$center,$counter,$width,10,$motif_color{$motif_id}) if($domains{$motif_id} or $arg{-d} =~ /all/); 
			$im->ellipse($margin+$centerx+$center,$counter,$width,10,$black) if($domains{$motif_id} or $arg{-d} =~ /all/);
			
		}
		elsif($shape =~ /p1/){
			pin1(\$im,$start,$end,$counter,$margin,$center,$motif_color{$motif_id});
		
		}		
		elsif($shape =~ /p2/){
			pin2(\$im,$start,$end,$counter,$margin,$center,$motif_color{$motif_id});
		
		}		
		elsif($shape =~ /p3/){
			pin3(\$im,$start,$end,$counter,$margin,$center,$motif_color{$motif_id});
		
		}		
		elsif($shape =~ /p4/){
			pin4(\$im,$start,$end,$counter,$margin,$center,$motif_color{$motif_id});
		
		}		
		elsif($shape =~ /p5/){
			pin5(\$im,$start,$end,$counter,$margin,$center,$motif_color{$motif_id});
		
		}		
		elsif($shape =~ /line/){
			line(\$im,$start,$end,$counter,$margin,$center,$motif_color{$motif_id});
		
		}		
		
		
	}
	
		
}

## Draw references

$counter += $step;

foreach $key (keys %motif_color){
	$counter += 25 if($domains{$key} or $arg{-d} =~ /all/);
	if($motif_shape{$key} =~/r/){ 
		$im->filledRectangle($margin+$adjust , $counter-5 ,191+$adjust , $counter+5 , $motif_color{$key} ) if($domains{$key} or $arg{-d} =~ /all/); 
		$im->rectangle($margin+$adjust , $counter-5 ,191+$adjust , $counter+5 , $black ) if($domains{$key} or $arg{-d} =~ /all/);
	}
	elsif($motif_shape{$key} =~/c/) { 
		$im->filledEllipse($margin+$adjust+5,$counter,15,10,$motif_color{$key}) if($domains{$key} or $arg{-d} =~ /all/); 
		$im->ellipse($margin+$adjust+5,$counter,15,10,$black) if($domains{$key} or $arg{-d} =~ /all/);
	}
	elsif($motif_shape{$key} =~/p1/){pin1(\$im, $adjust, 1+$adjust, $counter+8,$margin,0,$motif_color{$key})}
	elsif($motif_shape{$key} =~/p2/){pin2(\$im, $adjust, 1+$adjust, $counter-8,$margin,0,$motif_color{$key})}
	elsif($motif_shape{$key} =~/p3/){pin3(\$im, $adjust, 1+$adjust, $counter ,$margin,0,$motif_color{$key})}
	elsif($motif_shape{$key} =~/p4/){pin4(\$im, $adjust, 1+$adjust, $counter+5,$margin,0,$motif_color{$key})}
	elsif($motif_shape{$key} =~/p5/){pin5(\$im, $adjust, 1+$adjust, $counter-8,$margin,0,$motif_color{$key})}
	elsif($motif_shape{$key} =~/line/){line(\$im, $adjust, 1+$adjust, $counter,$margin,0,$motif_color{$key})}	
	$im->string(gdSmallFont,198+$adjust,$counter-5,"Motif $key ".$motif_seq{$key},$black) if($domains{$key} or $arg{-d} =~ /all/);
}

## make sure we are writing to a binary stream
binmode STDOUT;

## Convert the image to PNG and print it on standard output
print $im->png;	

sub pin1{
	my ($im,$start,$end,$counter,$margin,$center,$color) = @_;
	my $width = $end - $start;
	my $centerx = ($width/2)+$start;
	${$im}->filledEllipse($margin+$adjust+$centerx,$counter-15,7,7,$color); 
	${$im}->line($margin+$center+$centerx,$counter,$margin+$center+$centerx,$counter-10,$color);

}
sub pin2{
	my ($im,$start,$end,$counter,$margin,$center,$color) = @_;
	my $width = $end - $start;
	my $centerx = ($width/2)+$start;
	${$im}->filledEllipse($margin+$adjust+$centerx,$counter+15,7,7,$color); 
	${$im}->line($margin+$center+$centerx,$counter,$margin+$center+$centerx,$counter+10,$color);

}
sub pin3{
	my ($im,$start,$end,$counter,$margin,$center,$color) = @_;
	my $width = $end - $start;
	my $centerx = ($width/2)+$start;
	${$im}->filledEllipse($margin+$adjust+$centerx,$counter,12,12,$color); 
}
sub pin4{
	my ($im,$start,$end,$counter,$margin,$center,$color) = @_;
	my $width = $end - $start;
	my $centerx = ($width/2)+$start;
	${$im}->filledRectangle($margin+$center+$centerx-3,$counter-11,$margin+$center+$centerx+3,$counter-5,$color);
	${$im}->line($margin+$center+$centerx,$counter,$margin+$center+$centerx,$counter-5,$color);
 
}
sub pin5{
	my ($im,$start,$end,$counter,$margin,$center,$color) = @_;
	my $width = $end - $start;
	my $centerx = ($width/2)+$start;
	${$im}->filledRectangle($margin+$center+$centerx-3,$counter+5,$margin+$center+$centerx+3,$counter+11,$color);
	${$im}->line($margin+$center+$centerx,$counter,$margin+$center+$centerx,$counter+5,$color);
}
sub line{
	my ($im,$start,$end,$counter,$margin,$center,$color) = @_;
	my $width = $end - $start;
	my $centerx = ($width/2)+$start;
	${$im}->line($margin+$center+$centerx,$counter-3,$margin+$center+$centerx,$counter+3,$color);
}
