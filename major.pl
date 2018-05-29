#!/usr/bin/perl -w
use strict;

open(F1, "d:/code/pcrsim/genotypes.csv") or die("Could not open input file");

my $line;
my ($i,$j);
my @loci = ('D3','VWA','D16','D2','AMELO','D8','D21','D18','D19','TH','FGA');

for($j=0;$j<2;$j++){
	
	if($j==0){
		print "Major\n-----\n";
	}else{
		print "\nMinor\n-----\n";
	}
	my %target;
	my ($locus, $a1, $a2);
	
	for($i=0;$i<11;$i++){
		$line = <F1>;
		$line =~ s/[\r\n]+//g;
		($locus, $a1, $a2) = split(/,/,$line,3);
		if($a1!~/X/){
			if($a1>$a2){
				$target{$locus} = "$a2\,$a1";
			}else{
				$target{$locus} = "$a1\,$a2";
			}
		}else{
			$target{$locus} = "$a1\,$a2";
		}
	}
	
	foreach $locus (@loci){
		print "$locus,$target{$locus}\n";
	}
}
close(F1);
	
	