#!/usr/bin/perl -w
use strict;

my $line;
my $burnin = $ARGV[0];
my $sampleevery = $ARGV[1];
my $outputfile = $ARGV[2];

print "Burn In: $burnin\nSample every: $sampleevery\nOutputfile: $outputfile\n";

open(F1, "results.csv") or die("Couldn't open input file");
open(F2, ">$outputfile") or die("Couldn't open $outputfile for writing");

my $linecount = 0;
while($line=<F1> && $linecount<$burnin){
	$linecount++;
}

while($line=<F1>){
	if($linecount % $sampleevery == 0){
		print F2 $line;
	}
	$linecount++;
}

close(F1);
close(F2);
	
