#!/bin/bash
##
## procWTMM1d.sh


# NOTES: 
# i) The exponents tau(q) are computed for a list of moments q as defined 
# in the article (p.20). It would be probably simpler (for the discussion)
# to use a set of regularly spaced moments q, so that we don't have to
# focus on this question also. The commands below can be simply modified
# to fit the requirements of the article like in the first computation in 
# A.1 below (comment the commands between the flags EITHER and OR and 
# uncomment the commands between the flags OR and AND; modify similarly the 
# other cases).
# ii) The range of scale for the regression is r=[4,100], like in the
# article.
# iii) All the other parameters are "standard":
#   - The wavelet used for the analysis is the mexican hat (second 
#     derivative of wavelet)
#   - The range of scales of analysis is computed in regard to the size of
#     the input signal. The number of voices is fixed to 10.
# iv) We also compute the multifractal exponents tau(q)



DIRMF=..
DIRWTMM=../wtmm
RELDIRWTMM=../../wtmm
DIRBENCH=$DIRMF/benchmark
RELBINDIR=../../bin
WTMMPROCESSOR=wtmm_processor.exe
PROG=$RELBINDIR/$WTMMPROCESSOR



##  Ensemble A: 1x16384 (a single series with 16384 points)
SIZE=16386
N=1

cd $DIRBENCH/A
echo "Computing Distributions for signals of set A..."

## A.1 Log-Poisson h_\infty 1 = -0.5, D_\infty = 0
echo "         - A.1 Log-Poisson h_\infty = -0.5, D_\infty = 0"
# EITHER:
 # computation with the list of q as in the article: 
 echo "                     +++++++++++++++++++++++" > info.wtmm.setA
 echo " command: $WTMMPROCESSOR -f Log-Poisson.h-0.50-coD1.00_1D -N $N -regsc 4. 100. -tauq -dq 0." \
     >> info.wtmm.setA
 $PROG -f Log-Poisson.h-0.50-coD1.00_1D -N $N -regsc 4. 100. -tauq  -dq 0. \
     >> info.wtmm.setA
 mv Tq_wtmm* Dh_wtmm* $RELDIRWTMM/A/
# Note: we don't get significative difference
# OR:
# # standard computation with a list of moments q=[-8,8] regularly spaced (dq=0.2)
# echo "                     +++++++++++++++++++++++" > info.wtmm.setA
# echo " command: $WTMMPROCESSOR -f Log-Poisson.h-0.50-coD1.00_1D -N $N -regsc 4. 100. -q -8. 8. -tauq" \
#     >> info.wtmm.setA
# $PROG -f Log-Poisson.h-0.50-coD1.00_1D -N $N -regsc 4. 100.-q -8. 8. -tauq >> info.wtmm.setA
# mv Tq_wtmm* Dh_wtmm* $RELDIRWTMM/A/
# AND:
# computation through canonical approach: note unrelevant difference for such theoretical signal
echo >> info.wtmm.setA
echo "                    +++++++++++++++++++++++" >> info.wtmm.setA
echo " command: $WTMMPROCESSOR -f Log-Poisson.h-0.50-coD1.00_1D -N $N -regsc 4. 100. -dq 0. -tauq  -canon" \
    >> info.wtmm.setA
echo >> info.wtmm.setA
$PROG -f Log-Poisson.h-0.50-coD1.00_1D -N $N -regsc 4. 100. -tauq -canon \
    >> info.wtmm.setA
for f in Tq_wtmm* Dh_wtmm* 
do
  mv $f $RELDIRWTMM/A/can$f
done

## A.2 Log-Poisson h_\infty 1 = -0.25, D_\infty = 0.25
# computation with the list of q as in the article: 
echo "         - A.2 Log-Poisson h_\infty 1 = -0.25, D_\infty = 0.25"
echo >> info.wtmm.setA
 echo "                     +++++++++++++++++++++++" >> info.wtmm.setA
 echo " command: $WTMMPROCESSOR -f Log-Poisson.h-0.25-coD0.75_1D -N $N -regsc 4. 100. -tauq -dq 0." \
     >> info.wtmm.setA
 $PROG -f Log-Poisson.h-0.25-coD0.75_1D -N $N -regsc 4. 100. -tauq  -dq 0. \
     >> info.wtmm.setA
 mv Tq_wtmm* Dh_wtmm* $RELDIRWTMM/A/
# computation through canonical approach: note unrelevant difference for such theoretical signal
echo >> info.wtmm.setA
echo "                    +++++++++++++++++++++++" >> info.wtmm.setA
echo " command: $WTMMPROCESSOR -f Log-Poisson.h-0.25-coD0.75_1D -N $N -regsc 4. 100. -dq 0. -tauq  -canon" \
    >> info.wtmm.setA
echo >> info.wtmm.setA
$PROG -f Log-Poisson.h-0.25-coD0.75_1D -N $N -regsc 4. 100. -tauq -canon \
    >> info.wtmm.setA
for f in Tq_wtmm* Dh_wtmm* 
do
  mv $f $RELDIRWTMM/A/can$f
done

## A.3 Log-Normal h_m = 0.5, sigma_h = 1
# computation with the list of q as in the article: 
echo "         - A.3 Log-Normal h_m = 0.5, sigma_h = 1"
echo >> info.wtmm.setA
 echo "                     +++++++++++++++++++++++" >> info.wtmm.setA
 echo " command: $WTMMPROCESSOR -f Log-Normal.mean0.50-sigma1.00_1D -N $N -regsc 4. 100. -tauq -dq 0." \
     >> info.wtmm.setA
 $PROG -f Log-Normal.mean0.50-sigma1.00_1D -N $N -regsc 4. 100. -tauq  -dq 0. \
     >> info.wtmm.setA
 mv Tq_wtmm* Dh_wtmm* $RELDIRWTMM/A/
# computation through canonical approach: note unrelevant difference for such theoretical signal
echo >> info.wtmm.setA
echo "                    +++++++++++++++++++++++" >> info.wtmm.setA
echo " command: $WTMMPROCESSOR -f Log-Normal.mean0.50-sigma1.00_1D -N $N -regsc 4. 100. -dq 0. -tauq  -canon" \
    >> info.wtmm.setA
echo >> info.wtmm.setA
$PROG -f Log-Normal.mean0.50-sigma1.00_1D -N $N -regsc 4. 100. -tauq -canon \
    >> info.wtmm.setA
for f in Tq_wtmm* Dh_wtmm* 
do
  mv $f $RELDIRWTMM/A/can$f
done

## A.4 Log-Normal h_m = 0.25, sigma_h = 0.33
# computation with the list of q as in the article: 
echo "         - A.4 Log-Normal h_m = 0.25, sigma_h = 0.33"
echo >> info.wtmm.setA
 echo "                     +++++++++++++++++++++++" >> info.wtmm.setA
 echo " command: $WTMMPROCESSOR -f Log-Normal.mean0.25-sigma0.33_1D -N $N -regsc 4. 100. -tauq -dq 0." \
     >> info.wtmm.setA
 $PROG -f Log-Normal.mean0.25-sigma0.33_1D -N $N -regsc 4. 100. -tauq  -dq 0. \
     >> info.wtmm.setA
 mv Tq_wtmm* Dh_wtmm* $RELDIRWTMM/A/
# computation through canonical approach: note unrelevant difference for such theoretical signal
echo >> info.wtmm.setA
echo "                    +++++++++++++++++++++++" >> info.wtmm.setA
echo " command: $WTMMPROCESSOR -f Log-Normal.mean0.25-sigma0.33_1D -N $N -regsc 4. 100. -dq 0. -tauq  -canon" \
    >> info.wtmm.setA
echo >> info.wtmm.setA
$PROG -f Log-Normal.mean0.25-sigma0.33_1D -N $N -regsc 4. 100. -tauq -canon \
    >> info.wtmm.setA
for f in Tq_wtmm* Dh_wtmm* 
do
  mv $f $RELDIRWTMM/A/can$f
done

mv info.wtmm.setA $RELDIRWTMM/A/

################
##  Ensemble B: 10x4096 (10 series of 4096 points each)
SIZE=4096
N=10

cd ../B
echo "Computing Distributions for signals of set B..."

## B.1 Log-Poisson h_\infty = -0.5, D_\infty  = 0
# computation with the list of q as in the article: 
echo "         - B.1 Log-Poisson h_\infty = -0.5, D_\infty = 0"
echo >> info.wtmm.setB
 echo "                     +++++++++++++++++++++++" >> info.wtmm.setB
 echo " command: $WTMMPROCESSOR -f Log-Poisson.h-0.50-coD1.00_1D -N $N -regsc 4. 100. -tauq -dq 0." \
     >> info.wtmm.setB
 $PROG -f Log-Poisson.h-0.50-coD1.00_1D -N $N -regsc 4. 100. -tauq  -dq 0. \
     >> info.wtmm.setB
 mv Tq_wtmm* Dh_wtmm* $RELDIRWTMM/B
# computation through canonical approach: note unrelevant difference for such theoretical signal
echo >> info.wtmm.setB
echo "                    +++++++++++++++++++++++" >> info.wtmm.setB
echo " command: $WTMMPROCESSOR -f Log-Poisson.h-0.50-coD1.00_1D -N $N -regsc 4. 100. -dq 0. -tauq  -canon" \
    >> info.wtmm.setB
echo >> info.wtmm.setB
$PROG -f Log-Poisson.h-0.50-coD1.00_1D -N $N -regsc 4. 100. -tauq -canon \
    >> info.wtmm.setB
for f in Tq_wtmm* Dh_wtmm* 
do
  mv $f $RELDIRWTMM/B/can$f
done

## B.2 Log-Poisson h_\infty = -0.25, D_\infty  = 0.25
echo "         - B.2 Log-Poisson h_\infty 1 = -0.25, D_\infty = 0.25"
echo >> info.wtmm.setB
 echo "                     +++++++++++++++++++++++" >> info.wtmm.setB
 echo " command: $WTMMPROCESSOR -f Log-Poisson.h-0.25-coD0.75_1D -N $N -regsc 4. 100. -tauq -dq 0." \
     >> info.wtmm.setB
 $PROG -f Log-Poisson.h-0.25-coD0.75_1D -N $N -regsc 4. 100. -tauq  -dq 0. \
     >> info.wtmm.setB
 mv Tq_wtmm* Dh_wtmm* $RELDIRWTMM/B
# computation through canonical approach: note unrelevant difference for such theoretical signal
echo >> info.wtmm.setB
echo "                    +++++++++++++++++++++++" >> info.wtmm.setB
echo " command: $WTMMPROCESSOR -f Log-Poisson.h-0.25-coD0.75_1D -N $N -regsc 4. 100. -dq 0. -tauq  -canon" \
    >> info.wtmm.setB
echo >> info.wtmm.setB
$PROG -f Log-Poisson.h-0.25-coD0.75_1D -N $N -regsc 4. 100. -tauq -canon \
    >> info.wtmm.setB
for f in Tq_wtmm* Dh_wtmm* 
do
  mv $f $RELDIRWTMM/B/can$f
done

## B.3 Log-Normal h_m = 0.5, sigma_h = 1
echo "         - B.3 Log-Normal h_m = 0.5, sigma_h = 1"
echo >> info.wtmm.setB
 echo "                     +++++++++++++++++++++++" >> info.wtmm.setB
 echo " command: $WTMMPROCESSOR -f Log-Normal.mean0.50-sigma1.00_1D -N $N -regsc 4. 100. -tauq -dq 0." \
     >> info.wtmm.setB
 $PROG -f Log-Normal.mean0.50-sigma1.00_1D -N $N -regsc 4. 100. -tauq  -dq 0. \
     >> info.wtmm.setB
 mv Tq_wtmm* Dh_wtmm* $RELDIRWTMM/B
# computation through canonical approach: note unrelevant difference for such theoretical signal
echo >> info.wtmm.setB
echo "                    +++++++++++++++++++++++" >> info.wtmm.setB
echo " command: $WTMMPROCESSOR -f Log-Normal.mean0.50-sigma1.00_1D -N $N -regsc 4. 100. -dq 0. -tauq  -canon" \
    >> info.wtmm.setB
echo >> info.wtmm.setB
$PROG -f Log-Normal.mean0.50-sigma1.00_1D -N $N -regsc 4. 100. -tauq -canon \
    >> info.wtmm.setB
for f in Tq_wtmm* Dh_wtmm* 
do
  mv $f $RELDIRWTMM/B/can$f
done

## B.4 Log-Normal h_m = 0.25, sigma_h = 0.33
echo "         - B.4 Log-Normal h_m = 0.25, sigma_h = 0.33"
echo >> info.wtmm.setB
 echo "                     +++++++++++++++++++++++" >> info.wtmm.setB
 echo " command: $WTMMPROCESSOR -f Log-Normal.mean0.25-sigma0.33_1D -N $N -regsc 4. 100. -tauq -dq 0." \
     >> info.wtmm.setB
 $PROG -f Log-Normal.mean0.25-sigma0.33_1D -N $N -regsc 4. 100. -tauq  -dq 0. \
     >> info.wtmm.setB
 mv Tq_wtmm* Dh_wtmm* $RELDIRWTMM/B
# computation through canonical approach: note unrelevant difference for such theoretical signal
echo >> info.wtmm.setB
echo "                    +++++++++++++++++++++++" >> info.wtmm.setB
echo " command: $WTMMPROCESSOR -f Log-Normal.mean0.25-sigma0.33_1D -N $N -regsc 4. 100. -dq 0. -tauq  -canon" \
    >> info.wtmm.setB
echo >> info.wtmm.setB
$PROG -f Log-Normal.mean0.25-sigma0.33_1D -N $N -regsc 4. 100. -tauq -canon \
    >> info.wtmm.setB
for f in Tq_wtmm* Dh_wtmm* 
do
  mv $f $RELDIRWTMM/B/can$f
done

mv info.wtmm.setB $RELDIRWTMM/B/


 
################
##  Ensemble C: 100x1024 (100 series with 1024 points each)
SIZE=1024
N=100

cd ../C/
echo "Computing Distributions for signals of set C..."

## C.1 Log-Poisson h_\infty = -0.5, D_\infty = 0
echo "         - C.1 C.1 Log-Poisson h_\infty = -0.5, D_\infty = 0"
echo >> info.wtmm.setC
 echo "                     +++++++++++++++++++++++" >> info.wtmm.setC
 echo " command: $WTMMPROCESSOR -f Log-Poisson.h-0.50-coD1.00_1D -N $N -regsc 4. 100. -tauq -dq 0." \
     >> info.wtmm.setC
 $PROG -f Log-Poisson.h-0.50-coD1.00_1D -N $N -regsc 4. 100. -tauq  -dq 0. \
     >> info.wtmm.setC
 mv Tq_wtmm* Dh_wtmm* $RELDIRWTMM/C
# computation through canonical approach: note unrelevant difference for such theoretical signal
echo >> info.wtmm.setC
echo "                    +++++++++++++++++++++++" >> info.wtmm.setC
echo " command: $WTMMPROCESSOR -f Log-Poisson.h-0.50-coD1.00_1D -N $N -regsc 4. 100. -dq 0. -tauq  -canon" \
    >> info.wtmm.setC
echo >> info.wtmm.setC
$PROG -f Log-Poisson.h-0.50-coD1.00_1D -N $N -regsc 4. 100. -tauq -canon \
    >> info.wtmm.setC
for f in Tq_wtmm* Dh_wtmm* 
do
  mv $f $RELDIRWTMM/C/can$f
done

## C.2 Log-Poisson h_\infty = -0.25, D_\infty = 0.25
echo "         - C.2 Log-Poisson h_\infty 1 = -0.25, D_\infty = 0.25"
echo >> info.wtmm.setC
 echo "                     +++++++++++++++++++++++" >> info.wtmm.setC
 echo " command: $WTMMPROCESSOR -f Log-Poisson.h-0.25-coD0.75_1D -N $N -regsc 4. 100. -tauq -dq 0." \
     >> info.wtmm.setC
 $PROG -f Log-Poisson.h-0.25-coD0.75_1D -N $N -regsc 4. 100. -tauq  -dq 0. \
     >> info.wtmm.setC
 mv Tq_wtmm* Dh_wtmm* $RELDIRWTMM/C
# computation through canonical approach: note unrelevant difference for such theoretical signal
echo >> info.wtmm.setC
echo "                    +++++++++++++++++++++++" >> info.wtmm.setC
echo " command: $WTMMPROCESSOR -f Log-Poisson.h-0.25-coD0.75_1D -N $N -regsc 4. 100. -dq 0. -tauq  -canon" \
    >> info.wtmm.setC
echo >> info.wtmm.setC
$PROG -f Log-Poisson.h-0.25-coD0.75_1D -N $N -regsc 4. 100. -tauq -canon \
    >> info.wtmm.setC
for f in Tq_wtmm* Dh_wtmm* 
do
  mv $f $RELDIRWTMM/C/can$f
done

## C.3 Log-Normal h_m = 0.5, sigma_h = 1
echo "         - C.3 Log-Normal h_m = 0.5, sigma_h = 1"
echo >> info.wtmm.setC
 echo "                     +++++++++++++++++++++++" >> info.wtmm.setC
 echo " command: $WTMMPROCESSOR -f Log-Normal.mean0.50-sigma1.00_1D -N $N -regsc 4. 100. -tauq -dq 0." \
     >> info.wtmm.setC
 $PROG -f Log-Normal.mean0.50-sigma1.00_1D -N $N -regsc 4. 100. -tauq  -dq 0. \
     >> info.wtmm.setC
 mv Tq_wtmm* Dh_wtmm* $RELDIRWTMM/C
# computation through canonical approach: note unrelevant difference for such theoretical signal
echo >> info.wtmm.setC
echo "                    +++++++++++++++++++++++" >> info.wtmm.setC
echo " command: $WTMMPROCESSOR -f Log-Normal.mean0.50-sigma1.00_1D -N $N -regsc 4. 100. -dq 0. -tauq  -canon" \
    >> info.wtmm.setC
echo >> info.wtmm.setC
$PROG -f Log-Normal.mean0.50-sigma1.00_1D -N $N -regsc 4. 100. -tauq -canon \
    >> info.wtmm.setC
for f in Tq_wtmm* Dh_wtmm* 
do
  mv $f $RELDIRWTMM/C/can$f
done

## C.4 Log-Normal h_m = 0.25, sigma_h = 0.33
echo "         - C.4 Log-Normal h_m = 0.25, sigma_h = 0.33"
echo >> info.wtmm.setC
 echo "                     +++++++++++++++++++++++" >> info.wtmm.setC
 echo " command: $WTMMPROCESSOR -f Log-Normal.mean0.25-sigma0.33_1D -N $N -regsc 4. 100. -tauq -dq 0." \
     >> info.wtmm.setC
 $PROG -f Log-Normal.mean0.25-sigma0.33_1D -N $N -regsc 4. 100. -tauq  -dq 0. \
     >> info.wtmm.setC
 mv Tq_wtmm* Dh_wtmm* $RELDIRWTMM/C
# computation through canonical approach: note unrelevant difference for such theoretical signal
echo >> info.wtmm.setC
echo "                    +++++++++++++++++++++++" >> info.wtmm.setC
echo " command: $WTMMPROCESSOR -f Log-Normal.mean0.25-sigma0.33_1D -N $N -regsc 4. 100. -dq 0. -tauq  -canon" \
    >> info.wtmm.setC
echo >> info.wtmm.setC
$PROG -f Log-Normal.mean0.25-sigma0.33_1D -N $N -regsc 4. 100. -tauq -canon \
    >> info.wtmm.setC
for f in Tq_wtmm* Dh_wtmm* 
do
  mv $f $RELDIRWTMM/C/can$f
done

mv info.wtmm.setC $RELDIRWTMM/C/



################
# back...
cd $RELDIRWTMM/C/