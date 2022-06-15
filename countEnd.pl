#!/bin/perl

use strict;
use warnings;

unless(open(SAM, $ARGV[0])){die;}
unless(open(REF, $ARGV[1])){die;}
my $origSuf = '.primary.sam';
my $replSuf = '.primary.endcount.csv';
my $filename = $ARGV[0];
$filename =~ s/sam/endcount.csv/;
(my $logfilename = $ARGV[0]) =~ s/sam/endcount.log/;


unless(open(OUTPUT, ">$filename")){die;}
unless(open(LOG, ">$logfilename")){die;}

my @tRNAEndCount;
my $refSeqCount = 0;
while (<REF>){
	if ($_ =~ m/^>/) {
		chomp $_;
		my @splitRefHeader = split(' ', $_);
		(my $currentSeqName = $splitRefHeader[0]) =~ s/^\>//;
		my $currentSeqLength = $splitRefHeader[6];

		#initialize the tRNA end count table with 0, based on sequence length
		for (my $i = 0 ; $i < $currentSeqLength ; $i++){
			$tRNAEndCount[$refSeqCount][0] = $currentSeqName;
			$tRNAEndCount[$refSeqCount][1] = $currentSeqLength;
			$tRNAEndCount[$refSeqCount][2] = 0;
			$tRNAEndCount[$refSeqCount][$i+3] = 0;
		}

		#print "$currentSeqName\t$currentSeqLength\n";
		$refSeqCount++;
	}
}

#exit;

close REF;
#print "$refSeqCount\n"; exit;

#for (my $i = 0 ; $i < $refSeqCount ; $i++){
#	print "$tRNAEndCount[$i][0]\t$tRNAEndCount[$i][1]\n";
#}
#exit;

#my @SAMTable;
while (<SAM>){
	chomp $_;
	my @splitSAMLine = split('\t', $_);
	my $currentMappedRef = $splitSAMLine[2];
	my $currentMapStart = $splitSAMLine[3];
	my $currentCIGAR = $splitSAMLine[5];
	my $tempCIGAR = $currentCIGAR;
	(my $convNumCIGAR = $currentCIGAR) =~ s/[A-Z]/_/g;
	$convNumCIGAR =~ s/_$//;
	(my $convCodeCIGAR = $currentCIGAR) =~ s/[0-9]/_/g;
	$convCodeCIGAR =~ s/_+/_/g;
	$convCodeCIGAR =~ s/^_//;
	my @CIGARLengthSeq = split ('_', $convNumCIGAR);
	my @CIGARCodeSeq = split ('_', $convCodeCIGAR);

	#print "L:@CIGARLengthSeq C:@CIGARCodeSeq\t";

	my $CIGARCodeCount = @CIGARCodeSeq;
	my $CIGARsum = 0;

	if ($currentCIGAR =~ /N/) {
		$CIGARsum = -1;	
	} elsif ( @CIGARLengthSeq == @CIGARCodeSeq ) {
		for (my $i = 0; $i < @CIGARCodeSeq; $i++){
			if ($CIGARCodeSeq[$i] eq "S") {$CIGARsum += 0;}
			if ($CIGARCodeSeq[$i] eq "M") {$CIGARsum += $CIGARLengthSeq[$i];}
			if ($CIGARCodeSeq[$i] eq "D") {$CIGARsum += $CIGARLengthSeq[$i];}
			if ($CIGARCodeSeq[$i] eq "I") {$CIGARsum += 0;}
		}
	} else {
		$CIGARsum = -1;
	}
	
	my $endCoordinate = ($currentMapStart -1) + $CIGARsum;
	print LOG "$currentMappedRef\tS:$currentMapStart\t$currentCIGAR\t$convNumCIGAR\t$convCodeCIGAR\tCS:$CIGARsum\tE:$endCoordinate\tC:$CIGARCodeCount\n";
	#print "$currentMappedRef\tS:$currentMapStart\t$currentCIGAR\t$convNumCIGAR\t$convCodeCIGAR\tCS:$CIGARsum\tE:$endCoordinate\tC:$CIGARCodeCount\n";

	for (my $i = 0 ; $i < $refSeqCount ; $i++){
		if ($tRNAEndCount[$i][0] eq $currentMappedRef){
			#print "$i $tRNAEndCount[$i][0]\t$tRNAEndCount[$i][1]\n";
			if ( $endCoordinate > $tRNAEndCount[$i][1] ) { print "ERROR" ; die;}
			$tRNAEndCount[$i][2]++;
			$tRNAEndCount[$i][$endCoordinate+3]++;
		}
	}
	#print "\n";
}


print OUTPUT "Name,tRNALength,ReadCount,PostionalCount->\n";
for (my $i = 0 ; $i < $refSeqCount ; $i++){
	my $currentRefSeqLength = $tRNAEndCount[$i][1];
	print OUTPUT "$tRNAEndCount[$i][0]";
	for (my $j = 1 ; $j < $currentRefSeqLength+3 ; $j++){
		print OUTPUT ",$tRNAEndCount[$i][$j]";
	}
	print OUTPUT "\n";
}

close OUTPUT;
close LOG;
close SAM;
exit;
