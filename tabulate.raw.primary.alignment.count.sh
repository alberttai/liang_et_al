#!/bin/bash

cat ./ref/conv.hg19-mature-tRNAs.fa | grep \> | cut -d' ' -f1 > tRNAList.txt

filelist=$(ls *.primary.alignment.count.txt)
samplelist=${filelist//.tRNA.primary.alignment.count.txt/}
samplelist=$(echo $samplelist)
samplelist=${samplelist// /$'\t'}
header=$(echo -e "tRNA\t$samplelist")

umicollapsereadcount='rawReadCount'
for countfile in $(echo "$filelist") ; do
	samplereadcount=$(cat ${countfile/.tRNA.primary.alignment.count.txt/.process.log} | grep "unique sequence" | cut -d' ' -f8)
	umicollapsereadcount=$(echo -en "$umicollapsereadcount\t$samplereadcount")
done

counttable=$(echo -n "")

while read -r line ; do
	tRNA=${line#>}

	tRNAcount=''
	tRNAcountrow=$(echo -ne "${line#>}")
	for countfile in $(echo "$filelist") ; do
		sampleid=${countfile%.tRNA.primary.alignment.count.txt}
		#echo $sampleid $tRNA
		tRNAcount=$(cat $countfile | grep -P "\s$tRNA$")
		tRNAcount=$(echo $tRNAcount | cut -d' ' -f1)
		
		if [[ "$tRNAcount" == "" ]] ; then
			tRNAcountrow=$(echo -ne "$tRNAcountrow\t0")
		else
			tRNAcountrow=$(echo -ne "$tRNAcountrow\t$tRNAcount")
		fi
	done
	echo $line
	counttable=$(echo -e "$counttable\n$tRNAcountrow")

done < tRNAList.txt

echo "$header" > tRNARawReadCountRable.txt
echo "$umicollapsereadcount" >> tRNARawReadCountRable.txt
echo "$counttable" >> tRNARawReadCountRable.txt

cat tRNARawReadCountRable.txt
