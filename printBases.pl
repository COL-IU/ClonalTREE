#!/usr/bin/perl

my $pileupFile = $ARGV[0];

open (PILEUP, "$pileupFile") || die "cant load file";

while (<PILEUP>)
{
    chomp; #clear whitespaces
    @input = split("\t", $_); #splits each line with tabs
    $ref = $input[2]; 
    $lcref = lc($ref); #lowecase it
    $mut = $input[3];
    if (($mut =~ /\*/) or ($ref =~ /\*/))
    {
	next;
    }
    $input[4] =~ s/\Q.\E/$ref/g;   #$input[8] =~ s/\Q.\E/$ref/g;
    $input[4] =~ s/\Q,\E/$lcref/g; #$input[8] =~ s/\Q,\E/$lcref/g;
    $input[4] =~ s/\$//g;          #$input[8] =~ s/\$//g;
    #print "$input[3]";
    @bases = split("", $input[4]); #@bases = split("", $input[8]);
    
    #print @bases;
    #print "\n";
    $aCount=0;$ACount=0;
    $gCount=0;$GCount=0;
    $cCount=0;$CCount=0;
    $tCount=0;$TCount=0;
    $starCount=0;
    $indelCount=0;
    for($x=0;$x<scalar(@bases);$x++)
    {
	$indelSize = ();
	if ($bases[$x] =~ /\^/)
	{
	    $x++;
	    next;
	}
	#doesnt count indel
	if ($bases[$x] =~ /\+|-/)
	{
	    $x++;
	    $indelCount++;
	    while($bases[$x] =~ /(\d)/)
	    {
		$indelSize = $indelSize.$bases[$x];
		$x++;
	    }
	    for ($u=0;$u<$indelSize;$u++)
	    {
		$x++;
	    }
	}
	if ($bases[$x] =~ /A/){
	    $ACount++;
	}
	if ($bases[$x] =~ /a/){
	    $aCount++;
	}
	if ($bases[$x] =~ /C/){
	    $CCount++;
	}
	if ($bases[$x] =~ /c/){
	    $cCount++;
	}
	if ($bases[$x] =~ /G/){
	    $GCount++;
	}
	if ($bases[$x] =~ /g/){
	    $gCount++;
	}
	if ($bases[$x] =~ /T/){
	    $TCount++;
	}
	if ($bases[$x] =~ /t/){
	    $tCount++;
	}
	if ($bases[$x] =~ /\*/){
	    $starCount++;
	}
    }
    print join ("\t", $input[0], $input[1], $input[2], $input[3], $ACount, $aCount, $CCount, $cCount, $GCount, $gCount, $TCount, $tCount, $starCount, $indelCount), "\n";
    
}


