a="2.042"
b=`echo "2*$a" | bc -l`
for ix in `seq 0 $1`; do
  for iy in `seq 0 $2`; do
    for iz in `seq 0 $3`; do
      x=`echo "$ix*$b" | bc -l`
      y=`echo "$iy*$b" | bc -l`
      z=`echo "$iz*$b" | bc -l`
      echo "Li $x $y $z"
      x=`echo "$ix*$b+$a" | bc -l`
      y=`echo "$iy*$b+$a" | bc -l`
      z=`echo "$iz*$b" | bc -l`
      echo "Li $x $y $z"
      x=`echo "$ix*$b+$a" | bc -l`
      y=`echo "$iy*$b" | bc -l`
      z=`echo "$iz*$b+$a" | bc -l`
      echo "Li $x $y $z"
      x=`echo "$ix*$b" | bc -l`
      y=`echo "$iy*$b+$a" | bc -l`
      z=`echo "$iz*$b+$a" | bc -l`
      echo "Li $x $y $z"
      x=`echo "$ix*$b+$a" | bc -l`
      y=`echo "$iy*$b+$a" | bc -l`
      z=`echo "$iz*$b+$a" | bc -l`
      echo "H $x $y $z"
      x=`echo "$ix*$b+$a" | bc -l`
      y=`echo "$iy*$b" | bc -l`
      z=`echo "$iz*$b" | bc -l`
      echo "H $x $y $z"
      x=`echo "$ix*$b" | bc -l`
      y=`echo "$iy*$b+$a" | bc -l`
      z=`echo "$iz*$b" | bc -l`
      echo "H $x $y $z"
      x=`echo "$ix*$b" | bc -l`
      y=`echo "$iy*$b" | bc -l`
      z=`echo "$iz*$b+$a" | bc -l`
      echo "H $x $y $z"
    done  
  done
done
