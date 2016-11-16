# multifractal - Multifractal experiments
Tools for multifractal analysis of 1D (time-series) and 2D (images) signals

In the following, we describe the code source files that will enable 
you to reproduce the experiments of the article

 "Numerical methods for the estimation of singularity spectra on sampled 
 data: a comparative study",
 Turiel A., Pérez-Vicente C.J., and Grazzini, J.
 Journal of Computational Physics, 216(1):362-390, 2006.
 DOI: 10.1016/j.jcp.2005.12.004
 URL: http://www.sciencedirect.com/science/article/pii/S0021999105005565

You need to create your own hierarchy of directories, as follows:

          EXP -------- benchmark ------ A
                 |                 |
                 |                 |--- B
                 |                 |
                 |                 |--- C
                 |          
                 |          
                 |---- wtmm ----------- A     
                 |                 |
                 |                 |--- B
                 |                 |
                 |                 |--- C
                 |          
                 |          
                 |---- lastwave ------- A
                 |                 |
                 |                 |--- B
                 |                 |
                 |                 |--- C
                 |          
                 |          
                 |---- wavelab -------- A
                 |                 |
                 |                 |--- B
                 |                 |
                 |                 |--- C
                 |                 |
                 |                 |--- wavelib
                 |                 |
                 |                 |--- test
                 |          
                 |          
                 |---- fraclab -------- A
                                   |
                                   |--- B
                                   |
                                   |--- C
                                   |
                                   |--- fraclib
                                   |
                                   |--- test
 
Note that the name of the main directory (_EXP_ in the example above) 
doesn't matter whereas the others do.

benchmark
---------

The directory _benchmark/_ will contain all the data used in the various
experiments. 

In this directory, there is a script:

	buildbench.sh 
that enables to generate all the multifractal signals and to install 
them in the corresponding sub-directories _A/, B/ and C/_ once 
_multifractal_generator_ has been copied in _bin/_.
(note : in this file, change the name _multifractal_generator.exe_
into _multifractal_generator_)

wtmm
----

The directory _wtmm/_ will contain all the results of the WTMM analysis 
with our code.

In this directory, there is a script:

	procwtmm.sh
that enables to compute the singularity spectra of the signals in 
benchmark and to install the results of these estimations in the 
corresponding sub-directories _A/, B/ and C/_ once _wtmm_processor_ has 
been copied in bin/.
(note : in _procwtmm.sh_ 
   i/ first change the name _wtmm_processor.exe_ into _wtmm_processor_, 
   ii/ then read carefully the first lines of comments to have info about the processing, 
  iii/ the program needs approx. 20 mn to run with our data)

lastwave
--------

The directory _lastwave/_ stores the main program for analyzing the 
signals with the WTMM method implemented in the LastWave software.

In this directory, the file analysisLastWTMM1d.def contains the 
routines to process the analysis. To use the functions defined in 
this file, launch the software LastWave and type:

 	a> source analysisLastWTMM1d.def
("a>" is the prompt of lastwave). With the command:

	a> help AnalyzeSeriesLastWTMM 
you will have some information about the main function used to 
analyze our data.

If you want to analyze for instance the 10 signals Log-Poisson of 
size 4096 and with parameters h_infty=-0.5 and D_\infty=0, type the 
command:

	a> AnalyzeSeriesLastWTMM '../benchmark/B' 'B' \
		'Log-Poisson.h-0.50-coD1.00_1D' \
		10 1.5 8. 10. 'g2' {-200} -200 -200 
See the help of function _AnalyzeSeriesLastWTMM_ and the code for a 
complete description.

wavelab
-------

The directory _wavelab/_ contains the main program for analyzing the
signals with the WTMM method implemented in the WaveLab802 software.

  1. either you install the library WaveLab in your environment; then, 
  you  add the path of the directory WaveLab to the PATH variable of 
  Matlab, e.g. thanks to the command addpath (See the help of function 
  _procWaveWTMM1d_, below).
  2. you use the 'stand-alone' files provided in the sub-directory
  _wavelib/_ of the directory wavelab/. They are copies of the wavelab  
  necessary files. 

In Matlab, move to the directory _wavelab/_. Then, in order to perform 
the estimation for all the series, you should type:

	> procWaveWTMM1d('all','all');
(">" is the prompt of Matlab). The processing of all the series takes 
about 15mn. 
However, you may want to try first over small dataset to begin, so 
type for example:

	> procWaveWTMM1d('A',1);
to perform the analysis only for the multifractal serie stored in file 
_Log-Poisson.h-0.50-coD1.00_1D-N000.txt_ of set A. For more details, get 
the help of function _procWaveWTMM1d_ by typing:

	> help procWaveWTMM1d
The function performing the main computation is stored in the file 
_myAnalyzeSeriesWaveWTMM_; type

	> help myAnalyzeSeriesWaveWTMM
for help and have a look on the code inside this file.
Note that this method uses only the classical method to compute the 
Legendre transform (not the canonical one, like in lastwave).

You can have a look in the test file, whose purpose is to use different 
implementations (but all using the main WaveLab functions) for 
retrieving similar results.

fraclab
-------

The directory _fraclab/_ contains the main program for analyzing the 
signals with the WTMM method implemented in the FracLab software.

First of all, you should be sure that the fraclab functions are known 
in the Matlab environment. For this purpose, you have two solutions:

  1. either you install the library FracLab in your environment; then, 
  you need to add the path of the directory FracLab to the PATH variable 
  of Matlab, e.g. thanks to the command addpath (See the help of function 
  _procFracWTMM1d_, below).
  2. you use the 'stand-alone' files provided in the sub-directory 
  _fraclib/_ of the directory _fraclab/_. 

In Matlab, move to the directory _fraclab/_. Then, in order to perform 
the estimation for all the series, you should type:

	> procFracWTMM1d('all','all');
The processing of all the series takes about 25mn. 
However, you may want to try first over small dataset to begin, so type 
for example:

	> procFracWTMM1d('A',1);
to perform the analysis only for the multifractal stored in file 
_Log-Poisson.h-0.50-coD1.00_1D-N000.txt_ of set A. For more details, get 
the help of function _procFracWTMM1d_ by typing:

	> help procFracWTMM1d
The function performing the main computation is stored in the file 
_myAnalyzeSeriesFracWTMM_; type

	> help myAnalyzeSeriesFracWTMM
for help and have a look on the code.
Note that this method uses only the classical method to compute the 
Legendre transform.

We however recommand first to have a look on the _README_ of the 
sub-directory _test/_ and to perform the operations described therein. 
The aim of the described operations is to help you to understand how
FracLab works, to see the results it provides and eventually see some 
of its limitations (not the canonical one, like in lastwave).
