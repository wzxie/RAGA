#!/usr/bin/env perl 
use strict;
use warnings;
use Getopt::Long;
use SVG;

my $usage = <<USAGE;
Usage:
     $0 -f input.fa -o output.svg
	
OPTIONS
     -f                  # Input sequence file.
     -o                  # Set the output file name.
     -bar-color          # Set the bar color.    default: #1E90FF
     -split-length       # <int>    default: 100
     -canvas-size        # Set canvas size, height*width. default: 500*500
     -h/help             # Print the Usage.
	 
DESCRIPTION
     $0: The distribution and length of reads were analyzed.
USAGE
if (@ARGV==0){die $usage}
elsif($ARGV[0] =~ /-h/) {die $usage}

my ($inputfile1, $outputfile1, $splitLength, $bg_color, $rect_color, $width, $height, $left_gap, $right_gap, $down_gap, $up_gap, $angle, $canvas_size);
GetOptions(
    "f:s" => \$inputfile1,
    "o:s" => \$outputfile1,
	"split-length:i" => \$splitLength,
	"bar-color:s" =>\$rect_color,
	"canvas-size:s" =>\$canvas_size,
);
$splitLength ||= 100;
$inputfile1 ||= '';
$outputfile1 ||= '';
$bg_color ||="#F0F0F0";
$rect_color ||="#1E90FF";
$canvas_size ||="500*500";
{
      my @arr=split(/\*/,$canvas_size);
	  $height = $arr[0];
	  $width = $arr[1];
	  #print"$arr[0]\t$arr[1]\n";
}
$left_gap ||=50;
$right_gap ||=50;
$down_gap ||=50;
$up_gap ||=50;
$angle ||=30;
my $rotate="rotate($angle)";
my $conversion=3.14159265358979/180;
my $radians=$angle*$conversion;
my $font_size=9;
open IN, $inputfile1 or die "Can not open file $inputfile1, $!\n";
my (%seq, @seq_length, $seq_id);
while (<IN>) {
    chomp;
    if (m/^>(\S+)/) {
        $seq_id = $1;
    }
    else {
        $_ = uc($_);
        $seq{$seq_id} .= $_;
    }
}
close IN;
foreach my $id (keys %seq)
{
      my $seq = $seq{$id};
	  my $length=length($seq);
	  push(@seq_length,$length);
	  #if($length>59000){print"$id\n";}
}
my @sort_seql=sort{$b<=>$a}(@seq_length);
my $number1=int($sort_seql[0]/$splitLength);
if($sort_seql[0]%$splitLength!=0){$number1+=1;}

my %seq_length;
for(my $a1=1;$a1<=$number1;$a1++)
{
      my $a2=$a1*$splitLength;
	  $seq_length{$a2}=0;
}

foreach my $length (@seq_length)
{
      my $num1=int($length/$splitLength);
	  if($length%$splitLength!=0){$num1+=1;}
	  $num1*=$splitLength;
	  $seq_length{$num1}+=1;
}

my (@names, @scores);
foreach my $key1(sort{$a<=>$b}  keys %seq_length)
{
      my $key2=$key1-$splitLength;
	  my $key3="$key2".'~'."$key1";
	  my $len_num=$seq_length{$key1};
	  if($len_num!=0)
	  {
	      push(@names,$key3);
	      push(@scores,$len_num);
		  print"$key3\t$len_num\n";
	  }
}

my $svg= SVG->new( width=>$width, height=>$height);
my $h_ratio0=draw_rect();
draw_axis($h_ratio0);

open main1,">$outputfile1";
my $out1 = $svg->xmlify();
print main1 "$out1\n";
close main1;
print "\n--OK!--\n";


sub draw_rect{
      my $d_width=$width;
	  my $d_height=$height;
	  my $r_width=($width-$left_gap-$right_gap)/(scalar(@names)*1.5+0.5);
      my @sort_scores=sort{$b<=>$a}(@scores);
	  my $length1=length("$sort_scores[0]")-3;
	  my $powerofN=10**$length1;
	  if($powerofN<=1){$powerofN=5;}
	  my $max1=int($sort_scores[0]/$powerofN)*$powerofN;
	  if($sort_scores[0]%$powerofN!=0)
	  {
	         $max1+=$powerofN;
			 #print"$max1\n";
			 if($powerofN>5){
			 if($max1%(5*$powerofN)!=0){$max1+=5*$powerofN-$max1%(5*$powerofN);}}
	 }
	  #print"$sort_scores[0]\n";
	  #print"$length1\n$powerofN\n$max1\n";
      my $h_ratio=($d_height-$down_gap-$up_gap)/$max1;
	  my $rect_color0=$rect_color;
	  my $count=0;
	 $svg -> rectangle(
		  x        =>  $left_gap,
		  y        =>  $up_gap,
		  width =>  $d_width-$left_gap-$right_gap,
		  height=>  $d_height-$up_gap-$down_gap,
		  style   =>{
			  'fill'                  =>  "$bg_color",
			  'stroke'            =>  'black',
			  'stroke-width' =>  '0',
			  'fill-opacity'=>1, #0.5
		}
	);
	  foreach my $r_height (@scores)
	  {
	         $count+=0.5;
			 $svg -> rectangle(
	              x        =>  $left_gap+$r_width*$count,
	              y        =>  $d_height-$down_gap-$r_height*$h_ratio,
	              width =>  $r_width,
				  height=>  $r_height*$h_ratio,
	              style   =>{
				      'fill'                  =>  "$rect_color0",
					  'stroke'            =>  'black',
					  'stroke-width' =>  '0',
				}
	        );
			$count++;
	  }
	  $count=0;
	  foreach my $xbar (@names)
	  {
	         $count+=0.5;
			 my $fomt_gap=$font_size/2*length($xbar);
			 my $goal_x=$left_gap+$r_width*$count;
			 my $goal_y=$d_height-$down_gap+6;
			 my $rotate_x=abs($goal_x*cos($radians)+$goal_y*sin($radians));
			 my $rotate_y=abs($goal_y*cos($radians)-$goal_x*sin($radians));
			 $svg->text(
				 x=>$rotate_x, 
				 y=>$rotate_y,
				 transform=>$rotate,
				 style=>{
					 'font-family'=>"Calibri",#"Courier",
					 'stroke'=>'none',
					 'font-size'=>"$font_size",
				}
			)->cdata("$xbar");
			$count++;
	  }
	  return $h_ratio;
}

sub draw_axis{
     my $h_ratio=shift(@_);#$h_ratio=($d_height-$down_gap-$up_gap)/$max1;
	 $svg->line(
		 x1=>$left_gap,y1=>$up_gap,
		 x2=>$left_gap,y2=>$height-$down_gap,
		 style=>{
			 'stroke'=>"black",
			 'stroke-width'=>'0.5',
			 'stroke-opacity'=>0.8,   
		}
	);
	 $svg->line(
		 x1=>$left_gap,y1=>$height-$down_gap,
		 x2=>$width-$right_gap,y2=>$height-$down_gap,
		 style=>{
			 'stroke'=>"black",
			 'stroke-width'=>'0.5',
			 'stroke-opacity'=>0.8,   
		}
	);
	 my @sort_scores=sort{$b<=>$a}(@scores);
	 my $max1=$sort_scores[0];
	 my $length2=length($max1)-2;
	 my $ruler=10**$length2;
	 for(my $a2=$ruler;$a2<=$max1;$a2+=$ruler)
	 {
		 #print "$a2\n";
		 $svg->line(
			 x1=>$left_gap,y1=>$height-$down_gap-$a2*$h_ratio,
			 x2=>$width-$right_gap,y2=>$height-$down_gap-$a2*$h_ratio,
			 style=>{
				 'stroke'=>"#666666",
				 'stroke-width'=>'0.5',
				 'stroke-opacity'=>0.8,   
				 'stroke-dasharray'=>'5,5',
				 'stroke-dasharray-width'=>'1',
			}
		);
		 my $ybar = "$a2";
		 my $fomt_gap=$font_size/2*length($ybar);
		 $svg->text(
			 x=>$left_gap-$fomt_gap-2, y=>$height-$down_gap-$a2*$h_ratio+3,
			 style=>{
				 'font-family'=>"Calibri",#"Courier",
				 'stroke'=>'none',
				 'font-size'=>"$font_size",
			}
		)->cdata($ybar);
	 }
	 {
		 my $ybar = "0";
		 my $fomt_gap=$font_size/2*length($ybar);
		 $svg->text(
			 x=>$left_gap-$fomt_gap-2, y=>$height-$down_gap+3,
			 style=>{
				 'font-family'=>"Calibri",#"Courier",
				 'stroke'=>'none',
				 'font-size'=>"$font_size",
			}
		)->cdata($ybar);
	}
}
