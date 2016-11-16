++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The aim of the present document is to help you to understand how FracLab 
works, and to see the results it provides.

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Verify that the path to the fraclab toolbox has been added to the 
PATH variable of the environment.
This can be done thanks to the command ADDPATH, e.g. typing something 
like: 
   > addpath(genpath('MyPathToFraclab/Fraclab'));
where in my case MyPathToFraclab/Fraclab is the path of the directory 
containing the Fraclab functions.

Move to the directory test/ containing the test files
	> cd MyPathTocripts/test/
where MyPathTocripts is the path to the provided functions my#. 

This directory contains files which are renamed copies of the multifractal 
series generated thanks to the program multifractal_generator:
 - LogPoisson.A1.txt => Log-Poisson.h-0.50-coD1.00_1D-N000.txt,    set A
 - LogNormal.A3.txt  => Log-Normal.mean0.50-sigma1.00_1D-N000.txt, set A
 - LogPoisson.B1.txt => Log-Poisson.h-0.50-coD1.00_1D-N000.txt,    set B
 - LogPoisson.B2.txt => Log-Poisson.h-0.25-coD0.75_1D-N000.txt,    set B
 - LogNormal.B3.txt  => Log-Normal.mean0.50-sigma1.00_1D-N000.txt, set B
 - LogNormal.B4.txt  => Log-Normal.mean0.25-sigma0.33_1D-N000.txt, set B
 - LogPoisson.C1.txt => Log-Poisson.h-0.50-coD1.00_1D-N000.txt,    set C
 - LogNormal.C3.txt  => Log-Normal.mean0.50-sigma1.00_1D-N000.txt, set C
(These files have been renamed essentially because the FracLab toolbox 
doesn't support the file names with the character '-'), plus a file 
devilStaircase.txt storing the Cantor set.


++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			   TEST FRACLAB
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

NOTE: THE MAIN INTRUCTIONS TO FOLLOW ARE IN CAPITALS.

1/ Launch Fraclab
-----------------

Launch the fractal toolbox in a Matlab session with the COMMAND:
	> fltool;
(note: please don't forget the ";" mark).


1/ Load the signal
------------------

Select one of the signals in the current test/ directory by HITTING THE
BUTTON
	Load
in the main window. Then in the window appearing select the file and load
it in ASCII format.
You can eventually visuaize it thanks to the View button. 

Note:
i/ The fracLab toolbox doesn't support the file names with any '-' 
character.
ii/ to avoid errors signals, when choosing the menus of FL analysis, 
always verify that the signal you want to process is selected.


2/ Select the appropriate scheme of analysis
--------------------------------------------

First verify that the signal you want to process is selected, then 
CHOOSE THE MENU:	
	1D signals Multifractal Spectra
This menu allow you to compute the Legendre Spectrum and the Large 
Deviation Spectrum for 1D signals which are either functions or 
measures. In our cases, we are dealing with the estimation of Legendre
Spectrum (LS) over 1D signals, so, please, SELECT THE SUB-MENU:
	Functions
THEN:
 	Legendre Spectrum 
Then, there are again three algorithms for computing the LS: 
  - one based on the discrete wavelet transform (DWT), 
  - one on the continuous wavelet transform (CWT),
  - one on box counting. 
LAUNCH THE CWT BASED SUB-MENU: it corresponds exaclty to the WTMM  
method. A window called:
	"CWT Based Legendre Spectrum"
appears.


3/ Parameterize and launch the analysis
--------------------------------------

In order to have more control about the parameters of the analysis,
HIT THE BUTTON:
	Advanced compute. 
A window called:
	"CWT Based Legendre Spectrum Estimation"
appears.

Note: the signal which will be analyzed is the one highlighted when you 
press Advanced compute. Verify that it is the one appearing as an Input
Signal, otherwise change the variable in the main fraclab window and 
refresh the name in the current window.

In this window, you'll be able to set the main parameters pertaining to 
the computation, first, of the CWT:

  - fmin and fmax let you choose the minimum and maximum frequencies of 
    analysis. The default values are the ones yielding maximal span 
    compatible with the size of the signal.
    You may change the extreme frequencies either by typing values under 
    fmin and fmax, or by using the predefined values on the menus to 
    the right. 
    In the wtmm_processor, we used the variables MInScale and MaxScale, 
    which were also defined, by default, in regard to the size of the 
    signal. So no real difference in these parameters.
    LEAVE UNCHANGED AT FIRST TRY.

  - The Voices parameters governs the number of intermediate frequencies 
    at which the CWT coefficients will be computed. 
    LEAVE UNCHANGED AT FIRST TRY.
 
  - Checking the Mirror item will deal the border effects by mirroring 
    the signal at its extremities. Otherwise, zero-padding is used. 
    So UNSELECT IT AT FIRST TRY.

  - Choose the size and type of the wavelet: SELECT MEXICAN HAT

Once all the parameters that define the wavelet transform are chosen, HIT
THE BUTTON:
	Compute WT. 
The output signal is a matrix of size:
    "number of voices" x "size of the original signal". 
It should appear in the Input CWT box just below the Compute WT button.

In the lower part of the window, you may specify some new parameters in 
order to perform the WTMM: 

   - Qmin, Qmax, and # of Q's: the function T(q), from whih LS is computed  
     as the Legendre transform  is estimated for q ranging in [Qmin, Qmax], 
     and # of Q's values are computed in this interval. 
     INCREASE THE RANGE OF Q's AND THEIR # (51 VALUES IN [-5,5] FOR 
     INSTANCE)

   - the Regression Range: the range of scales over which the tau(q) are 
     estimated. By default, the computation can be made for a interactive
     choice of this range.
     If you uncheck Specify Regression Range, the label Full Regression 
     Range appears instead, and the estimate will be performed on the whole 
     range of scales as defined by the number of Voices. 
     LEAVE THE bUTTON CHECKED TO TEST SEVERAL REGRESSIONS MANUALLY; THIS 
     WAY, PERFORM SEVERAL ESTIMATIONS.

   - the Type of Regression from the usual choice. 
     LEAVE UNCHANGED (LEAST SQUARES)

   - the Maxima Box: this performs the maxima selection in WTMM method. It 
     selects, at each scale, instead of a mean value, the largest 
     coefficients in given neighbourhoods.
     LEAVE ALWAYS THIS BOX MARKED.


4/ Play to produce different estimations
----------------------------------------

In the case you leaved Regression Range checked, a graphic window will 
appear after you HIT COMPUTE: You'll see a bunch of curves, one curve per 
value of q. 
For each q, T(q) is estimated as the slope of the the best linear fit of 
the corresponding curve. Thus, you will want to select, with the cross
shaped cursor, a range of scales (in abscissa), where an approximately 
linear behaviour holds for all or most of the curves. 

Once you have selected such a range, two new graphs will appear to the 
right: 
  - on the lower one, you'll see the estimate of T(q), 
  - on the upper one will be the Legendre spectrum D(h). 

Press return on the keyboard to finish. 

!!!
As they say in the documentation: 
	"You may experiment with other ranges until you are satisfied  
         with the result"
That means...? You decide when the results are good? not the algorithm?
!!!


++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			TEST "MY#FUNCTIONS#"
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


If necessary add the directory of the functions MYCWTPART and MYCWTSPEC, 
MYCONTWTPART and MYCONTWTSPEC to path variable of Matlab:
	> addpath '../'
In the directory test/, TYPE: 
	> help testFracFunctions
to test some of these functions.

For example:
	> testFracFunctions(1);
will enable to process the file Log-Poisson.h-0.50-coD1.00_1D-N000.txt of 
the dataset A.
	> testFracFunctions(5);
will enable to process the file Log-Normal.mean0.50-sigma1.00_1D-N000.txt 
of the dataset B.

Compare the results with those obtained thanks to the FracLab toolbox.




