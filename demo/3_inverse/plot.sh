#!/bin/bash

set -e

p=../../plugin

${p}/alltone.py -0 init-361.dat -P 26 -b 250 output/allinone.txt -a avrgModel.txt -Sk

${p}/senker.py -W 200.0 -f f.txt -c c.txt -w  16:287,120:720 avrgModel.txt -Sk -E c=3100
${p}/senker.py -W 200.0 -f f.txt -c c.txt -w 161:510,350:700 avrgModel.txt -Sk -E c=4000
