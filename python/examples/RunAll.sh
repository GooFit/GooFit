#!/usr/bin/env bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

set -evx
cd $DIR
python -V
time python exponential.py
time python plot_2d.py
time python addition.py
time python convolution.py
time python DP4.py
time python product.py
time python TDDP4.py
time python SigGen.py
time python dalitz.py
time python chisquare.py
time python simpleFit.py
time python zachFit.py --reduce 10
#time python pipipi0.py
set +evx
