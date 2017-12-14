#!/bin/sh

# Group carriers into buckets
# Usage
#
# BucketFile file BucketSize
# 
# The input file syntax should be:
# e -8.53033e-16 400.0 280.0 2.0e-9
#
# Geneva 3 Nov 2017

fnm=$1
Bsz=$2

rm -f BUCKETED
Nlines=`wc -l $fnm | awk '{print $1}'`
Nbkt=`echo "round($Nlines/$Bsz)" | calc -p`
for i in `seq 1 $Nbkt`
do

   ist=`echo "$Bsz*($i-1)+1" | calc -p`
   iend=`echo "$ist+$Bsz-1" | calc -p`
   awk -v i1=$ist -v i2=$iend -v bsz=$Bsz '{\
        if ( NR>=i1 && NR<=i2 ) {
          qsum+=$2;
	  ax+=$3;
	  ay+=$4;
        } 
	ax=ax/bsz;
	ay=ay/bsz;
	print($1," ",qsum," ",ax," ",ay,$5);
   }' $fnm >> BUCKETED
  read a
done
