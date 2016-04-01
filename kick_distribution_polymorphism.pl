#!/usr/bin/perl

#*--<<Definision>>-----------------------------------------*
#PGID         Visualize_alignment_gene.pl
#Kind of PG   Main PG
#Create Date  2010/07/28
#
#	驕ｺ莨晏ｭ仙腰菴阪〒alignment繧恥ng繧､繝｡繝ｼ繧ｸ縺ｫ縺吶ｋ縲�#
#		
#   Comandline  gff file
#               multiple fasta file
#				output image file
#               pileup file
#*---------------------------------------------------------*
#***********************************************************
# use
#************************************************************
use strict;
use warnings;
use CGI;
use Template;
use FindBin;
#***********************************************************
# variable declaration
#***********************************************************
my $reference;
my $pileup;
my $base = 10;
my $output;
my $directory_len;
my $directory;
my @list;
my $script_dir = $FindBin::Bin;
#***********************************************************
# Main Coading
#***********************************************************
$reference            = $ARGV[0];
$pileup               = $ARGV[1];
$output               = $ARGV[2];

#make a directory set alignment file
$directory_len = length($output);
$directory = substr($output,0,$directory_len-4);
$directory = $directory."_files";
system("mkdir $directory");

system("perl $script_dir/ForDistriSNP.pl '$reference' '$pileup' '$base' '$directory'");
system("perl $script_dir/makeImage_imlib_distriSNP_cairo.pl '$reference' '$base' '$directory'");

&data;
&html;
sub data{

	my $dh;
	my $file;
	my $i=0;
	opendir $dh,$directory or die "Directory for results is invalid.";
	my @files = sort { $a cmp $b } readdir($dh);
	while ( my $file = shift @files ) {
#	while (my $file = readdir $dh) {
		if($file =~ /\.png$/){
			my $explain     = substr($file,0,index($file,"\."));
			$list[$i] = ({png =>"$file", explain => $explain},);
			$i++;
		}
	}
	
}
sub html{

	my $cgi = CGI->new;
	my $title = 'Distribution of DNA polymorphisms in chromosomes';
	my $direction = 'This number of DNA polymorphisms in each chromosome is shown in brackets. The x-axis represents the physical distance along each chromosome, split into 500-kb windows. The orange lines represent regions in which no DNA polymorphisms were detected. The y-axis indicates the common logarithm of number of DNA polymorphisms.';

	system("cp $script_dir/style.css $directory");

	my $output_tmp = '';
	my $template = Template->new({
   		INCLUDE_PATH => $script_dir,
	});
	$template->process(
  	'templete_distributionSNP.html',
  	{
    	style => $directory."/"."style.css",
    	title => $title, 
		direction => $direction,
    	list => \@list,
  	},
  	\$output_tmp,
	);

	open(OUTPUT, ">$output");
	print $cgi->header();
	print OUTPUT $output_tmp;
	close OUTPUT;

}	
