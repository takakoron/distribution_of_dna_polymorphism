#! /usr/bin/perl

#*--<<Definision>>-----------------------------------------*
#PGID         makeCoverByPileUp_100kb_each.pl
#Kind of PG   Main PG
#Create Date  2010/07/28
#
#	Distribution of DNApolymorphism....preprocess 
#		
#   Comandline  reference multiple fasta file(chromosome name‚Íchr01‚ðŽg—p)
#               pileup file
#				output text file
#*---------------------------------------------------------*
#***********************************************************
# use
#***********************************************************
use strict;
#***********************************************************
# constant
#***********************************************************
my $spLeng = "500000";     ##REWRITABLE   split length
my $shift  = "500000";     ##REWRITABLE    shift length
#***********************************************************
# variable
#***********************************************************
my @chrNo;
my %chrLeng;
my $line;
my $reference;
my $pileup;
my $base;
my $output_dir;
#***********************************************************
# Main Coading
#***********************************************************
$reference            = $ARGV[0];
$pileup               = $ARGV[1];
$base                 = $ARGV[2];
$output_dir           = $ARGV[3];

open (FASTA,    "<$reference")      or die "$reference : No such file or directory\n";
my $i=0;
my $seq;
while($line = <FASTA>){
	if($line =~ /^>/){
		chomp $line;
		if($seq ne ""){
			$chrLeng{$chrNo[$i]} = length($seq);
			$seq = "";
			$i++;
		}
		my $chr_tmp;
		if(index($line," ") > 0){
			$chrNo[$i] = substr($line,4,index($line," ")-4);
			$chr_tmp = substr($line,1,index($line," ")-1);
		}else{
			$chrNo[$i] = substr($line,4);
			$chr_tmp = substr($line,1);
		}
		
		my $output_tmp = $output_dir.'/'.$chr_tmp.'.pileup';
		system("grep $chr_tmp $pileup > $output_tmp ")
	}else{
		$seq .= $line;
	}
}
$chrLeng{$chrNo[$i]} = length($seq);
$seq = "";
			
foreach my $i(0..$#chrNo) {
	my %sp = ();
	my $limit = "";
	if($chrLeng{$chrNo[$i]}%$shift > 0) {
		$limit = int($chrLeng{$chrNo[$i]}/$shift);
	}else {
		$limit = int($chrLeng{$chrNo[$i]}/$shift) - 1;
	}
		
	foreach my $j(0..$limit) {
		$sp{$j} = 0;
	}

	my $filename = $output_dir.'/chr'.$chrNo[$i].'.pileup';
	open(FILE,$filename) || die "could not open $filename\n";
	while(<FILE>) {
		my @line = split('\t',$_);
		if(exists($sp{int($line[1]/$shift)})) {
			$sp{int($line[1]/$shift)}++; 
		} 
 		if(int($line[1]/$shift) > 0){   ## more than 0 as 1st...
			if($line[1] >= (int($line[1]/$shift) -1)*$shift + $shift && $line[1] <= (int($line[1]/$shift) -1)*$shift + $spLeng) {
				$sp{int($line[1]/$shift) - 1}++;
			}
		}
	}
	close(FILE);
	my @keylist = keys %sp;
	@keylist = sort {$a <=> $b} @keylist;

	my $filename = $output_dir.'/chr'.$chrNo[$i].'.txt';
	open(OUT,">$filename");
	foreach my $j(0..$#keylist) {
		my $no = $j + 1;
		my $to = $j*$shift + $spLeng;

		if($j == $#keylist) {
			if($keylist[$j]*$shift+$spLeng > $chrLeng{$chrNo[$i]}) {
				$to = $chrLeng{$chrNo[$i]};
			}
		}
		if($sp{$j}  == 0) {
			if($base eq ""){
				print OUT $no."\t".($j*$shift + 1)."\t".$to."\t".$sp{$j}."\t".log($sp{$j} + 1)."\n"; 
			}else{
				print OUT $no."\t".($j*$shift + 1)."\t".$to."\t".$sp{$j}."\t".log($sp{$j} + 1)/log($base)."\n"; 
			}
		}else {
			if($base eq ""){
				print OUT $no."\t".($j*$shift + 1)."\t".$to."\t".$sp{$j}."\t".(int(log($sp{$j})*100+0.5)/100)."\n";
			}else{		
				print OUT $no."\t".($j*$shift + 1)."\t".$to."\t".$sp{$j}."\t".(int(log($sp{$j})/log($base)*100+0.5)/100)."\n"; 
			}
		}

		#print OUT $no."\t".($j*$shift + 1)."\t".$to."\t".$sp{$j}."\n";
	}
	close(OUT);
}
