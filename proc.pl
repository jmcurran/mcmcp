open(F1, "results10k.csv") or die("Couldn\'t open results file");

while($line=<F1>){
    $line=~s/[\r\n]+//g;
    @tokens = split(/[ ]+/,$line,44);
    @majorgenotypes = @tokens[0..12];
    @minorgenotypes = @tokens[13..25];
    $keymajor = join(' ',@majorgenotypes);
    $countsmajor{$keymajor}++;
    $keyminor = join(' ',@minorgenotypes);
    $countsminor{$keyminor}++;
}
close(F1);

open(F2,">outmajor.csv") or die("Couldn\'t open the outfile ");
foreach $key (sort keys %countsmajor){
    print F2 "$key\,$countsmajor{$key}\n";
}

close(F2);

open(F2,">outminor.csv") or die("Couldn\'t open the outfile ");
foreach $key (sort keys %countsminor){
    print F2 "$key\,$countsminor{$key}\n";
}

close(F2);
