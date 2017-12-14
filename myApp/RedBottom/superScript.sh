#! /bin/bash

for i in `seq -w 1 4`
do
  touch `date +%Y%m%d`-$i.txt
  ../MfgTRACSFit $i 2015-09-01_15-21-40_PIN_7859_1.8_4_W2_J6-1_Vbias100-500.txt.root Config.TRACS_10._1. >> `date +%Y%m%d`-$i.txt
  rm file*
  mv output.root output-$i.root
done
