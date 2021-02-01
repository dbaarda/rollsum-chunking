#!/usr/bin/bash
dir=${1:-.}
[[ -d $dir ]] || mkdir $dir
for alg in weibull0 weibull1 weibull2 fastcdc; do
  ./chunker.py $alg $dir >$dir/$alg.txt
done
./graphs.py $dir
