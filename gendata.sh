#!/usr/bin/bash
CHUNKERS="chunker weibull1 weibull2 weibullt1 weibullt2 nc1 nc2 nc3 rc4"
dir=${1:-.}
[[ -d $dir ]] || mkdir $dir
for alg in $CHUNKERS; do
  ./chunker.py $alg $dir >$dir/$alg.txt
done
./graphs.py $dir
