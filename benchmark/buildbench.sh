#!/bin/bash
##
## buildbench.sh
DIRMF=../msm-wtmm
DIRBENCH=.
# $HOME/Multifractal/MFestimation
# MSMDIR=$DIRMF/msm-wtmm
BINDIR=$DIRMF/bin
MSMPROCESSOR=multifractal_generator.exe
# msm_processor.exe
PROG=$BINDIR/$MSMPROCESSOR


DSPACE=1

################
##  Ensemble A: 1x16384 (a single series with 16384 points)
SIZE=16386
N=1

TYPE=0
## A.1 Log-Poisson h_\infty 1 = -0.5, D_\infty = 0
# default choice for h_\infty and D_\infty 
echo "            +++++++++++++++++++++++" > info.setA
echo "  command: $MSMPROCESSOR -dim $SIZE -N $N -type $TYPE -ver " >> info.setA
echo >> info.setA
$PROG -dim $SIZE -N $N -type $TYPE -ver >> info.setA
#
## A.2 Log-Poisson h_\infty 1 = -0.25, D_\infty = 0.25
HINF=-0.25
CODINF=0.75
echo "            +++++++++++++++++++++++" >> info.setA
echo "  command: $MSMPROCESSOR -dim $SIZE -N $N -type $TYPE -hinf $HINF -Codinf $CODINF -ver " \
    >> info.setA
echo >> info.setA
$PROG -dim $SIZE -N $N -type $TYPE -ver  -hinf $HINF -Codinf $CODINF >> info.setA

TYPE=1
## A.3 Log-Normal h_m = 0.5, sigma_h = 1
# default choice for h_m and sigma_h
echo "            +++++++++++++++++++++++" >> info.setA
echo "  command: $MSMPROCESSOR -dim $SIZE -N $N -type $TYPE -ver " >> info.setA
echo >> info.setA
$PROG -dim $SIZE -N $N -type $TYPE -ver >> info.setA
#
## A.4 Log-Normal h_m = 0.25, sigma_h = 0.33
HMEAN=0.25
SIGMA=0.33
echo "            +++++++++++++++++++++++" >> info.setA
echo "  command: $MSMPROCESSOR -dim $SIZE -N $N -type $TYPE -mean $HMEAN -sigma $SIGMA -ver " \
    >> info.setA
echo >> info.setA
$PROG -dim $SIZE -N $N -type $TYPE -mean $HMEAN -sigma $SIGMA -ver >> info.setA

mv info.setA $DIRBENCH/A
mv Log*dat $DIRBENCH/A
mv Log*hdr $DIRBENCH/A
mv Log*txt $DIRBENCH/A
mv Dh* $DIRBENCH/A


################
##  Ensemble B: 10x4096 (10 series of 4096 points each)
SIZE=4096
N=10
#
TYPE=0
## B.1 Log-Poisson h_\infty = -0.5, D_\infty  = 0
# default choice for h_\infty and D_\infty 
echo "            +++++++++++++++++++++++" > info.setB
echo "  command: $MSMPROCESSOR -dim $SIZE -N $N -type $TYPE -ver " >> info.setB
echo >> info.setB
$PROG -dim $SIZE -N $N -type $TYPE -ver >> info.setB
#
## B.2 Log-Poisson h_\infty = -0.25, D_\infty  = 0.25
HINF=-0.25
CODINF=0.75
echo "            +++++++++++++++++++++++" >> info.setB
echo "  command: $MSMPROCESSOR -dim $SIZE -N $N -type $TYPE -hinf $HINF -Codinf $CODINF -ver " \
    >> info.setB
echo >> info.setB
$PROG -dim $SIZE -N $N -type $TYPE -ver  -hinf $HINF -Codinf $CODINF >> info.setB

TYPE=1
## B.3 Log-Normal h_m = 0.5, sigma_h = 1
# default choice for h_m and sigma_h
echo "            +++++++++++++++++++++++" >> info.setB
echo "  command: $MSMPROCESSOR -dim $SIZE -N $N -type $TYPE -ver " >> info.setB
echo >> info.setB
$PROG -dim $SIZE -N $N -type $TYPE -ver >> info.setB
#
## B.4 Log-Normal h_m = 0.25, sigma_h = 0.33
HMEAN=0.25
SIGMA=0.33
echo "            +++++++++++++++++++++++" >> info.setB
echo "  command: $MSMPROCESSOR -dim $SIZE -N $N -type $TYPE -mean $HMEAN -sigma $SIGMA -ver " \
    >> info.setB
echo >> info.setB
$PROG -dim $SIZE -N $N -type $TYPE -mean $HMEAN -sigma $SIGMA -ver >> info.setB

mv info.setB $DIRBENCH/B
mv Log*dat $DIRBENCH/B
mv Log*hdr $DIRBENCH/B
mv Log*txt $DIRBENCH/B
mv Dh* $DIRBENCH/B


 
################
##  Ensemble C: 100x1024 (100 series with 1024 points each)
SIZE=1024
N=100
#
TYPE=0
## C.1 Log-Poisson h_\infty = -0.5, D_\infty = 0
# default choice for h_\infty and D_\infty 
echo "            +++++++++++++++++++++++" > info.setC
echo "  command: $MSMPROCESSOR -dim $SIZE -N $N -type $TYPE -ver " >> info.setC
echo >> info.setC
$PROG -dim $SIZE -N $N -type $TYPE -ver >> info.setC

## C.2 Log-Poisson h_\infty = -0.25, D_\infty = 0.25
HINF=-0.25
CODINF=0.75
echo "            +++++++++++++++++++++++" >> info.setC
echo "  command: $MSMPROCESSOR -dim $SIZE -N $N -type $TYPE -hinf $HINF -Codinf $CODINF -ver " \
    >> info.setC
echo >> info.setC
$PROG -dim $SIZE -N $N -type $TYPE -ver  -hinf $HINF -Codinf $CODINF >> info.setC


TYPE=1
## C.3 Log-Normal h_m = 0.5, sigma_h = 1
# default choice for h_m and sigma_h
echo "            +++++++++++++++++++++++" >> info.setC
echo "  command: $MSMPROCESSOR -dim $SIZE -N $N -type $TYPE -ver " >> info.setC
echo >> info.setC
$PROG -dim $SIZE -N $N -type $TYPE -ver >> info.setC
#
## C.4 Log-Normal h_m = 0.25, sigma_h = 0.33
HMEAN=0.25
SIGMA=0.33
echo "            +++++++++++++++++++++++" >> info.setC
echo "  command: $MSMPROCESSOR -dim $SIZE -N $N -type $TYPE -mean $HMEAN -sigma $SIGMA -ver " \
    >> info.setC
echo >> info.setC
$PROG -dim $SIZE -N $N -type $TYPE -mean $HMEAN -sigma $SIGMA -ver >> info.setC

mv info.setC $DIRBENCH/C
mv Log*dat $DIRBENCH/C
mv Log*hdr $DIRBENCH/C
mv Log*txt $DIRBENCH/C
mv Dh* $DIRBENCH/C
