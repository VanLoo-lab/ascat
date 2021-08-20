#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;
use File::Basename;
use POSIX;

my %options=();
getopts("s:p:c:", \%options);

print "Please provide chromosome sizes\n" and die if !defined $options{s};
print "Please provide SNP pos file for array\n" and die if !defined $options{p};
if(!defined $options{c}){
    $options{c}=1;
    print "No number of cores defined. Default=1\n";
}

#store chr sizes
my %chrLength=();

open(IN, "<",$options{s}) or die "Couldn't find chromosome sizes file ".$options{s};
while(<IN>){
    chomp($_);
    my @tab=split(/\t+/, $_);
    $chrLength{$tab[0]}=$tab[1];
}
close(IN);


my @windows = (12,25,50,100,250,500,1000,2500,5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000);

open(IN, "<",$options{p}) or die "Couldn't find SNP pos platform file ".$options{p};
my $probes = @{[<IN>]};
close(IN);

my $lines=$probes*@windows;

my $counter=0;

open(OUT, ">", basename($options{p}, ",txt")."_".$counter.".bed");
open(IN, "<",$options{p}) or die "Couldn't find SNP pos platform file ".$options{p};
while(<IN>){
    chomp($_);
    if($_ !~ m/^\t/gi){
    my @tab=split(/\t+/, $_);
    
    foreach my $w(@windows){
        my $start=$tab[2]-$w;
        my $stop=$tab[2]+$w;
        
        #adjust for chromosome borders
        if($start<0){
            $start=0;
        }
        if($stop>($chrLength{"chr".$tab[1]})){
            $stop=$chrLength{"chr".$tab[1]};
        }
        
        if($w==12){
            print OUT "chr".$tab[1]."\t".$start."\t".$stop."\t".$tab[0]."\t".$tab[2]."\t".($w*2+1)."\n";
        }
        else{
           print OUT "chr".$tab[1]."\t".$start."\t".$stop."\t".$tab[0]."\t".$tab[2]."\t".($w*2)."\n";
        }
        $counter++;
        
        if($counter % ceil($lines/$options{c}) == 0){
            close(OUT);
            open(OUT, ">", basename($options{p}, ",txt")."_".int(($counter/($lines/$options{c}))+0.5).".bed");
        }
    }
    }
}
close(IN);
close(OUT);
