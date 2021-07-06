#!/bin/bash
err=.000001

for q in 5
do
  for l in 46
  do
    echo $l
    # python3.8 matgen.py $l 20 $(($q**2)) "$l"x20_q"$q"_c.txt
    for i in {4..5}
    do
      echo $i
      run sebd.exe $l 20 $q 1 $err 1 "$l"x20_q"$q"_c.txt o"$l"x20_q"$q"_e-6_"$i".txt 1 1
    done
  done
done
for q in 6
do
  for l in {22..46..8}
  do
    echo $l
    python3.8 matgen.py $l 20 $(($q**2)) "$l"x20_q"$q"_c.txt
    for i in {1..5}
    do
      echo $i
      run sebd.exe $l 20 $q 1 $err 1 "$l"x20_q"$q"_c.txt o"$l"x20_q"$q"_e-6_"$i".txt 1 1
    done
  done
done
echo "done"
