#!/bin/bash

set -e

p=../../plugin

${p}/alltone.py -0 ../data/init-361.dat -P 26 -b 250 output/allinone.txt \
  -a avrgModel.txt -Sk

${p}/spec.py +b data +p real -W 200.0 \
  -f ../data/f.txt -c ../data/c.txt -d ../data/361snSpec.txt -t \
  -0 ../data/init-361.dat -2 avrgModel.txt \
  -w 16:287,120:720 161:510,350:700 -S

${p}/senker.py -W 200.0 -f ../data/f.txt -c ../data/c.txt \
  -w  16:287,120:720 avrgModel.txt -Sk -E c=3100
${p}/senker.py -W 200.0 -f ../data/f.txt -c ../data/c.txt \
  -w 161:510,350:700 avrgModel.txt -Sk -E c=4000
