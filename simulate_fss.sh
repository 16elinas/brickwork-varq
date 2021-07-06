#!/bin/bash
for q in {7..8}
do
  for l in {24..32..2}
  do
      err=.000001
    # echo $err
    run sebd.exe $l 20 $q 1 $err 1 "$l"x20_q"$q"_c.txt o"$l"x20_q"$q"_e-6.txt 1 1
  done
done

for l in {18..24..2}
do
  err=.000001
  # echo $err
  run sebd.exe $l 20 8 1 $err 1 "$l"x20_q8_c.txt o"$l"x20_q8_e-6.txt 1 1
done


echo "done"

# run sebd.exe 22 20 7 1 .000001 1 22x20_q7_c.txt o22x20_q7_e-6.txt 1 1
