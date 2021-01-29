#!/usr/bin/bash

for alg in weibull0 weibull1 weibull2 fastcdc; do
  ./chunker.py $alg >data/$alg.txt
done
