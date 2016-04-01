#! /usr/bin/perl -w

#*--<<Definision>>-----------------------------------------*
#PGID         perl makeImage.pl
#Kind of PG   Main PG
#Create Date  2010/07/28
#
#	Distribution of DNApolymorphism....preprocess 
#		
#   Comandline  reference multiple fasta file
#               pileup file
#				output text file
#*---------------------------------------------------------*
#***********************************************************
# use
#***********************************************************
use strict;
use Data::Dumper;
use Cairo;
use Gtk2;
#***********************************************************
# constant
#***********************************************************
my @chrNo;
my %chrLeng;

#***********************************************************
# variable
#***********************************************************
my $reference;
my $base;
my $output_dir;
#***********************************************************
# Main Coading
#***********************************************************
$reference            = $ARGV[0];
$base                 = $ARGV[1];
$output_dir           = $ARGV[2];

open (FASTA,    "<$reference")      or die "$reference : No such file or directory\n";
my $i=0;
my $seq="";
while(my $line = <FASTA>){
	if($line =~ /^>/){
		chomp $line;
		if($seq ne ""){
			$chrLeng{$chrNo[$i]} = length($seq);
			$seq = "";
			$i++;
		}
		if(index($line," ") > 0){
			$chrNo[$i] = substr($line,4,index($line," ")-4);
		}else{
			$chrNo[$i] = substr($line,4);
		}
	}else{
		$seq .= $line;
	}
}
$chrLeng{$chrNo[$i]} = length($seq);
$seq = "";

foreach my $h(0..$#chrNo) {

	my $file = $output_dir.'/chr'.$chrNo[$h].'.txt';
	my $outFig = $output_dir."/chr".$chrNo[$h].".png";
	my %snp = ();
	my %code = ();

	undef my @order;
	my $winLengDefo = 0;
	open(FILE,$file) || die "could not open $file\n";
	while(<FILE>) {
		if($_ !~ /^\#\#/) {
			chomp;
			my @line = split('\t',$_);    ## window pos, cover rate, # of snp
			push @order,($line[1] - 1);   ##convert to strating position from 0 to 1
			if($#order  == 0) {
				$winLengDefo = $line[2] - $line[1] + 1; 
			}
			$snp{($line[1] - 1)} = $line[4];   ## log number of snps
 			$code{($line[1] -1)} = $line[3];   ## real number of snps
 			## mod to snps/Kbp 20090522
			#$snpNo{$line[0]} = $line[2];   ##FIX later
		}
	}
	close(FILE);


	my $marginL = 26;
	#my $marginR = 80;
	my $marginR = 50;
	my $marginA = 6;
	my $marginB = 18;

	my $maxCov = 4;   ## mod   ##height
	my $Cov = 20;       ## mod

	my $horiz = scalar keys %snp;
	my $horizW = 7;   ## horizonal width

	my $x = $marginL + $horiz*$horizW + $marginR;
	my $y = $marginA + $maxCov*$Cov + $marginB;

###############################################################################################
	my $surface = Cairo::ImageSurface->create ('argb32',$x,$y);
	my $cr = Cairo::Context->create ($surface);
	$cr->rectangle(1,1,$x,$y);
	$cr->set_source_rgb (1,1,1);
	$cr->fill;
	my $layout = Gtk2::Pango::Cairo::create_layout( $cr );
###############################################################################################

##title
	my $x1 = 0;
	my $y1 = 0;

##vartical line
	$x1 = $marginL;
	my $x2 = $marginL;
	$y1 = $marginA + $maxCov*$Cov;
	my $y2 = $marginA; 

###############################################################################################
   $cr->set_source_rgb(0,0,0);
   $cr->move_to ($x1,$y1);   #put start point to draw
   $cr->line_to ($x2,$y2);
   $cr->stroke;
###############################################################################################
   
	$x1 = $marginL -12;  ##mod   scale text for vartical line
	foreach my $i(0..$maxCov) {
		if($i % 1 == 0) {
			$y1 = $marginA + ($maxCov - $i)*$Cov -8;
###############################################################################################			
			$cr->move_to ($x1,$y1);
			my $context = $layout->get_context;			
			$cr->set_font_size(9);
			my $layout_text = "<span font_family='Courier' color='black'>$i</span>";
			$layout->set_markup($layout_text);
			Gtk2::Pango::Cairo::show_layout($cr, $layout);
###############################################################################################				
		}
	}

##horizonal line
	$x1 = $marginL;
	$y1 = $marginA + $maxCov*$Cov;
	#$x2 = $marginL + $horiz*$horizW + 1;
	#$y2 = $marginA + $maxCov*$Cov;
	$x2 = ($marginL + $horiz*$horizW + 1) - $marginL;
	$y2 = 3;  ##mod 20090609
###############################################################################################	
	$cr->set_source_rgb(0,0,225);      # blue
	$cr->rectangle($x1,$y1,$x2,$y2);   #put start point to draw
	$cr->fill;
###############################################################################################
   
	my $scale = 10000000;   ##scale of horizonal (10Mb)
	my $scaleNo = 1;
	foreach my $i(0..$#order) {    ##draw SNPs
		my @line = split('\t',$snp{$order[$i]});
		if($snp{$order[$i]} != -1) {
			my $snp =  int($snp{$order[$i]}*$Cov + 0.5);
			if($snp{$order[$i]} > $maxCov){    ##4 hight
				$snp = $maxCov*$Cov;   ##mod 20091105
			}
			$x1 = $marginL + $horizW*$i + 1;
			$x2 = $marginL + $horizW*$i + 1;
			
			if($i < $#order) {
				$x2 = $marginL + $horizW*($i+1) - ($marginL + $horizW*$i + 1);  ## MUST MOD 4 last colamn
			}else{   ## $i == $#order
				if($order[$i] + $winLengDefo >= $chrLeng{$chrNo[$h]}) {
					$x2 = int((($chrLeng{$chrNo[$h]} - $order[$i])/$winLengDefo*$horizW)+ 0.5);
				}
			}


			# $y1 = $marginA + $maxCov*$Cov - 1;
			# $y2 = $marginA + $maxCov*$Cov - $snp;
			$y1 = $marginA + $maxCov*$Cov - $snp; 
			$y2 = ($marginA + $maxCov*$Cov) - ($marginA + $maxCov*$Cov - $snp);
###############################################################################################
			$cr->set_source_rgb(0,0,205);      # blue
			$cr->rectangle($x1,$y1,$x2,$y2);   #put start point to draw
			$cr->fill;
###############################################################################################
		}

		if($code{$order[$i]} == 0) {    ##real number of SNPs
			$x1 = $marginL + $horizW*$i + 1;
			$x2 = $marginL + $horizW*($i+1) - ($marginL + $horizW*$i + 1);
			$y1 = $marginA + $maxCov*$Cov;
			$y2 = 3;
###############################################################################################	
			$cr->set_source_rgb(255,0,0);      # red
 			$cr->rectangle($x1,$y1,$x2,$y2);   #put start point to draw
			$cr->fill;
###############################################################################################
		}

		if($order[$i] + 1 > $scale*$scaleNo) {
			$x1 = $marginL + $horizW*($i-1); 
			$y1 = $marginA + $maxCov*$Cov + 1;
			my $text = $scale*$scaleNo/1000000;
###############################################################################################
			$cr->move_to ($x1,$y1);
			my $context = $layout->get_context;			
			my $layout_text = "<span font_family='Courier New' color='black'>$text</span>";
			$layout->set_markup($layout_text);
			Gtk2::Pango::Cairo::show_layout($cr, $layout);			
###############################################################################################
			$scaleNo++;
		}
		if($i == $#order) {
			$x1 = $marginL + $horiz*$horizW + 8;
			# $y1 = $marginA + $maxCov*$Cov -12;
			$y1 = $marginA + $maxCov*$Cov -1;
			my $text = "\(Mb\)";
			
###############################################################################################
			$cr->move_to ($x1,$y1);
			my $context = $layout->get_context;			
			$cr->set_font_size(11);
			my $layout_text = "<span font_family='Courier New' color='black' size='10000'>$text</span>";
			$layout->set_markup($layout_text);
			Gtk2::Pango::Cairo::show_layout($cr, $layout);			
###############################################################################################
		}
	}

	undef @order;


	binmode STDOUT;
###############################################################################################			
	$x1 = -4;
	$y1 = $marginA + $Cov - 15;
	$cr->move_to (0,80);
	my $context = $layout->get_context;			
	$cr->rotate(4.7);  #これでうまく縦書きっつーか横書きにできるのかな～？？？
	my $layout_text = "<span font_family='Courier' color='black' size='10000'>log(SNPs)</span>";
	$layout->set_markup($layout_text);
	Gtk2::Pango::Cairo::show_layout($cr, $layout);
###############################################################################################		

	$cr->show_page;
	$surface->write_to_png ($outFig);
	
#	undef $image;

}
