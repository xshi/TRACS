bnm=`basename $1 .txt`
grep "NumberDs in last Z step with Height" $1 | awk '{print $8,"",$9*100./36872}' | tr -d ':' | sort > $bnm.dat
