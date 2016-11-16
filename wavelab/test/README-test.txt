++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The aim of the present document is to help you to illustrate the use of 
some WaveLab functions and to compare with different approaches derived 
from it. 

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Move to the directory test/ containing the test files
	> cd MyPathTocripts/test/
where MyPathTocripts is the path to the provided functions my#. 

Verify that the path to the wavelab toolbox has been added to the 
PATH variable of the environment.
This can be done thanks to the command ADDPATH, e.g. typing something 
like: 
   > addpath(genpath('MyPathToWavelab/Wavelab'));
where in my case MyPathToWavelab/Wavelab is the path of the directory 
containing the Wavelab functions, or by adding the pass to the sub-dir
wavelib/ (upper directory):
   > addpath(genpath('../wavelib'));


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


++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			TEST "MY#FUNCTIONS#"
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


In the directory test/, TYPE: 
	> help testWaveFunctions
to test some of these functions.

For example:
	> testWaveFunctions(1);
will enable to process the file Log-Poisson.h-0.50-coD1.00_1D-N000.txt of 
the dataset A.
	> testWaveFunctions(5);
will enable to process the file Log-Normal.mean0.50-sigma1.00_1D-N000.txt 
of the dataset B.




