#!/usr/bin/env perl 
use strict;
use warnings;
use Getopt::Long;
use SVG;

my $usage = <<USAGE;
Usage:
     $0 -f input.bed -l genome_length.txt -o output.svg
	
OPTIONS
     -f                  # Input bed file.
     -l                  # Input genome-length file.
     -o                  # Set the output file name.
     -line-color         # set the line color  default: #1E90FF
     -split-length       # <int>    default: 100
     -canvas-size        # Set canvas size, height*width. default: 500*500
     -h/help             # Print the Usage.
	 
DESCRIPTION
     $0: The distribution and length of reads were analyzed.
USAGE
if (@ARGV==0){die $usage}
elsif($ARGV[0] =~ /-h/) {die $usage}

my ($inputfile1, $inputfile2, $outputfile1, $splitLength, $bg_color, $rect_color, $width, $height, $left_gap, $right_gap, $down_gap, $up_gap, $angle, $font_size, $canvas_size);
GetOptions(
               "f:s" => \$inputfile2,
               "l:s" => \$inputfile1,
               "o:s" => \$outputfile1,
    "split-length:i" => \$splitLength,
      "line-color:s" => \$rect_color,
     "canvas-size:s" => \$canvas_size,
);
$splitLength ||= 100;
$inputfile1 ||= '';
$outputfile1 ||= '';
$bg_color ||="#F0F0F0";
$rect_color ||="#1E90FF";
$canvas_size ||="500*500";
{
      my @arr = split(/\*/,$canvas_size);
      $height = $arr[0];
       $width = $arr[1];
}
$left_gap ||=80;
$right_gap ||=0;
$down_gap ||=0;
$up_gap ||=50;
$angle ||=30;
$font_size ||=15;
my $rotate="rotate($angle)";
my $conversion=3.14159265358979/180;
my $radians=$angle*$conversion;

my (%seq,@seqlength,%site);
my $contig_num=0;
open IN1, $inputfile1 or die "Can not open file $inputfile1, $!\n";
while (<IN1>) {
     chomp;
	 my @arr=split /\t/;
	 $seq{$arr[0]}=$arr[2];
         my $NA_l='NA';
         $site{$arr[0]}->{$NA_l}=$NA_l;
	 push(@seqlength,$arr[2]);
	 $contig_num++;
}
close IN1;
my @sort_seql=sort{$b<=>$a}(@seqlength);
my $maxlength=$sort_seql[0];

open IN2, $inputfile2 or die "Can not open file $inputfile2, $!\n";
while (<IN2>) {
     chomp;
	 my @arr=split/\t/;
	 $site{$arr[0]}->{$arr[1]}=$arr[2];
}
close IN2;


my $svg= SVG->new( width=>$width, height=>$height);
draw_genome();
draw_line();

open main1,">$outputfile1";
my $out1 = $svg->xmlify();
print main1 "$out1\n";
close main1;
print "\n--OK!--\n";


sub draw_genome{
      my $max_width=$width-$left_gap-2;
	  my $g_height=($height-$up_gap-2)/($contig_num*2);
	  my $ratio=$max_width/$maxlength;
	  my $count=0;
	  $svg -> rectangle(
		  x        =>  $left_gap,
		  y        =>  $up_gap,
		  width =>  $max_width,
		  height=>  $height-$up_gap-1,
		  style =>{
			  'fill'         =>  "$bg_color",
			  'stroke'       =>  'black',
			  'stroke-width' =>  '0',
			  'fill-opacity'=>1, #0.5
		}
	 );
	 $svg->line(
		 x1=>$left_gap,y1=>$up_gap,
		 x2=>$width-2,y2=>$up_gap,
		 style=>{
			 'stroke'=>"#666666",
			 'stroke-width'=>'1',
			 'stroke-opacity'=>1,
		}
	);
		 #Draw coordinate axes
		 my $t_postfix="Mb";
		 my $valueSpan= 5000000; 
		 my $scaleXX;
		 my $scaleYY;			
		 if(($maxlength)>=1000000000){$t_postfix="Gb";}
		 elsif(($maxlength)>=1000000){$t_postfix="Mb";}
		 elsif(($maxlength)>=1000){$t_postfix="Kb";}
		 elsif(($maxlength)>=1){$t_postfix="bp";}
		 if(($maxlength)>=1000000000){$valueSpan=1000000000;$scaleYY=1;$scaleXX=$valueSpan/10;}
		 elsif(($maxlength)>=100000000){$valueSpan=50000000;$scaleYY=50;$scaleXX=$valueSpan/5;}
		 elsif(($maxlength)>=10000000){$valueSpan=5000000;$scaleYY=5;$scaleXX=$valueSpan/5;}
		 elsif(($maxlength)>=1000000){$valueSpan=1000000;$scaleYY=1;$scaleXX=$valueSpan/10;}
		 elsif(($maxlength)>=100000){$valueSpan=50000;$scaleYY=50;$scaleXX=$valueSpan/5;}
		 elsif(($maxlength)>=10000){$valueSpan=5000;$scaleYY=5;$scaleXX=$valueSpan/5;}
		 elsif(($maxlength)>=1000){$valueSpan=1000;$scaleYY=1;$scaleXX=$valueSpan/10;}
		 elsif(($maxlength)>=100){$valueSpan=50;$scaleYY=50;$scaleXX=$valueSpan/5;}
		 elsif(($maxlength)>=10){$valueSpan=5;$scaleYY=5;$scaleXX=1;}			 
		 elsif(($maxlength)>=1){$valueSpan=1;$scaleYY=1;$scaleXX=1;}
		 my $t = 0;
		 for (my $i=0;$i<=$maxlength;$i+=$scaleXX)
		 {
			 my $i_pXX=$i*$ratio;
			 if(($i)%($valueSpan) ==0){
				 my $xbar = "$t"."$t_postfix";
				 $svg->text(
					 x=>$left_gap +$i_pXX, y=>$up_gap -5,
					 style=>{
						 'font-family'=>"Calibri",#"Courier",
						 'stroke'=>'none',
						 'font-size'=>$font_size-2,
					}
				)->cdata($xbar);
				 $t+=$scaleYY;
				 $svg->line(
					 x1=>$i_pXX+$left_gap , y1=>$up_gap -5,
					 x2=>$i_pXX+$left_gap , y2=>$up_gap,
					 style=>{
						 'stroke'=>'black',
						 'stroke-width'=>'0.5',
					}
				);
				 $svg->line(
					 x1=>$i_pXX+$left_gap , y1=>$up_gap,
					 x2=>$i_pXX+$left_gap , y2=>$height-1,
					 style=>{
						 'stroke' => 'black',
						 'stroke-dasharray'=>'3,3',
						 'stroke-width'=>'0.5',
					}
				);
			}
			 else
			 {
				 $svg->line(
					 x1=>$i_pXX+$left_gap , y1=>$up_gap -3,
					 x2=>$i_pXX+$left_gap , y2=>$up_gap,
					 style=>{
						 'stroke'=>'black',
						 'stroke-width'=>'0.5',
						}
				);	
				 $svg->line(
					 x1=>$i_pXX+$left_gap , y1=>$up_gap,
					 x2=>$i_pXX+$left_gap , y2=>$height-1,
					 style=>{
						 'stroke' => '#cccccc',
						 'stroke-dasharray'=>'3,3',
						 'stroke-width'=>'0.5',
					}
				);
			}
		}
	 foreach my $key1(sort keys %seq)
	 {
		  my $key2=$seq{$key1};
		  print"$key1\t$key2\n";
		  $count++;
		  $svg -> rectangle(
			  	    x        =>  $left_gap,
				    y        =>  $up_gap+$count*$g_height,
				    rx       =>  $g_height/2,
				    ry       =>  $g_height/2,
				    width    =>  $key2*$ratio,
				    height   =>  $g_height,
				    style    =>{
						 'fill'         =>  "#FFFFFF",
						 'stroke'       =>  'black',
						 'stroke-width' =>  '0.5',
					}
				);
		  my $s_length=length($key1);
		  $svg->text(
		  	     x    => $left_gap-$s_length*$font_size/2-5, 
			     y    => $up_gap+$count*$g_height+$g_height/2+$font_size/3,
			     style=>{
			   	     'font-family'=>"Calibri",#"Courier",
				     'stroke'     =>'none',
				     'font-size'  =>"$font_size",
				    }
		            )->cdata("$key1");
		  $count++;
	}
}

sub draw_line{
      my $max_width=$width-$left_gap-2;
	  my $g_height=($height-$up_gap-2)/($contig_num*2);
	  my $ratio=$max_width/$maxlength;
	  my $count=0;
	  foreach my $key1(sort keys %seq)
	 {
		  my $key2=$site{$key1};
                  #print"$key1\n";
		  $count++;
		  foreach my $key3(keys %{$key2})
		  {
				 my $key4=$site{$key1}->{$key3};
                                 if($key4 ne 'NA')
                                 {
			                $svg -> rectangle(
						  x     =>  $left_gap+$key3*$ratio,
						  y     =>  $up_gap+$count*$g_height+1,
						  width => ($key4-$key3)*$ratio,#*10,
						  height=>  $g_height-2,
						  style =>{
							  'fill' => "$rect_color",#"black",#
							  'stroke' => 'black',
							  'stroke-width' =>  '0',
							  'fill-opacity'=>1,
						}
					);
                                 }
		  }
		  $count++;
	}
         $count=0;
         foreach my $key1(sort keys %seq)
         {
                  my $key2=$seq{$key1};
                  $count++;
                  $svg -> rectangle(
                                    x        =>  $left_gap,
                                    y        =>  $up_gap+$count*$g_height,
                                    rx       =>  $g_height/2,
                                    ry       =>  $g_height/2,
                                    width    =>  $key2*$ratio,
                                    height   =>  $g_height,
                                    style    =>{
                                                 'fill'         =>  'none',
                                                 'stroke'       =>  'black',
                                                 'stroke-width' =>  '0.5',
                                        }
                                );
                  $count++;
        }
}
