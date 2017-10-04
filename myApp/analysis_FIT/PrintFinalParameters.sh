cat > TOLERANCE << EOF
10.
1.
0.1
0.01
EOF

cat > WEIGHT << EOF
1.
1.e-1
1.e-2
1.e-3
1.e-4
1.e-5
1.e-6
1.e-7
EOF

rm -f REPORT
for TOLE in `cat TOLERANCE`
do
  
  for WGT in `cat WEIGHT`
  do
  
  echo "Tolerance=" $TOLE "  WEIGHT=" $WGT      >> REPORT
  ~/bin/Greplog log_"$TOLE"_"$WGT" | tail -n 15 | grep -v "# ext." >> REPORT
  echo "======================================================================" >> REPORT
  
  done
  
  
done
