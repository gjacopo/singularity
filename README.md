singularity
===========

Numerical methods for multifractal singularity analysis of 1D (time-series) signals
---

**Description**

The code source files will enable you to reproduce the experiments on multifractal singularity analysis detailed in [Turiel _et al._'s article](References) on _numerical methods for the estimation of singularity spectra on sampled data: a comparative study_ referenced below (_you can cite this source code or the reference's doi: [10.1016/j.jcp.2005.12.004](http://dx.doi.org/10.1016/j.jcp.2005.12.004)_).

**Usage**

*settings*

You need to create your own hierarchy of directories, as follows:

          EXP -------- benchmark ------ A
                 |                 |--- B
                 |                 |--- C   
                 |          
                 |---- wtmm ----------- A     
                 |                 |--- B
                 |                 |--- C        
                 |          
                 |---- lastwave ------- A
                 |                 |--- B
                 |                 |--- C
                 |          
                 |---- wavelab -------- A
                 |                 |--- B
                 |                 |--- C
                 |                 |--- wavelib
                 |                 |--- test
                 |          
                 |---- fraclab -------- A
                                   |--- B
                                   |--- C
                                   |--- fraclib
                                   |--- test
 
Note that the name of the main directory (_EXP_ in the example above) doesn't matter whereas the others do.

*benchmark*

The directory `benchmark/` will contain all the data used in the various experiments. 

In this directory, there is a script:

	buildbench.sh 
that enables to generate all the multifractal signals and to install them in the corresponding sub-directories `A/`, `B/` and `C/` once `multifractal_generator` has been copied in `bin/` (note : in this file, change the name `multifractal_generator.exe` into `multifractal_generator`).

*wtmm*

The directory `wtmm/` will contain all the results of the WTMM analysis with our code.

In this directory, there is a script:

	procwtmm.sh
that enables to compute the singularity spectra of the signals in  benchmark and to install the results of these estimations in the corresponding sub-directories `A/`, `B/` and `C/` once `wtmm_processor` has been copied in `bin/`.
(note : in `procwtmm.sh`, 
   i/ first change the name `wtmm_processor.exe` into `wtmm_processor`, 
   ii/ then read carefully the first lines of comments to have info about the processing, 
  iii/ the program needs approx. 20 mn to run with our data)

*lastwave*

The directory `lastwave/` stores the main program for analyzing the signals with the WTMM method implemented in the [_LastWave_](LastWave) software.

In this directory, the file `analysisLastWTMM1d.def` contains the routines to process the analysis. To use the functions defined in this file, launch the software _LastWave_ and type:

 	a> source analysisLastWTMM1d.def
("a>" is the prompt of lastwave). With the command:

	a> help AnalyzeSeriesLastWTMM 
you will have some information about the main function used to analyze our data.

If you want to analyze for instance the 10 signals Log-Poisson of  size 4096 and with parameters `h_infty=-0.5` and `D_\infty=0`, type the command:

	a> AnalyzeSeriesLastWTMM '../benchmark/B' 'B' 'Log-Poisson.h-0.50-coD1.00_1D' \
		10 1.5 8. 10. 'g2' {-200} -200 -200 
See the help of function `AnalyzeSeriesLastWTMM` and the code for a complete description.

*wavelab*

The directory `wavelab/` contains the main program for analyzing the signals with the WTMM method implemented in the [_WaveLab802_][WaveLab850] software.

  1. either you install the library _WaveLab_ in your environment; then, you  add the path of the directory _WaveLab_ to the `PATH` variable of _Matlab_, _e.g._ thanks to the command addpath (See the help of function 
  `procWaveWTMM1d`, below).
  2. you use the 'stand-alone' files provided in the sub-directory
  `wavelib/` of the directory `wavelab/`. They are copies of the _Wavelab_  
  necessary files. 

In _Matlab_, move to the directory `wavelab/`. Then, in order to perform the estimation for all the series, you should type:

	> procWaveWTMM1d('all','all');
(">" is the prompt of Matlab). The processing of all the series takes about 15mn. However, you may want to try first over small dataset to begin, so type for example:

	> procWaveWTMM1d('A',1);
to perform the analysis only for the multifractal serie stored in file `Log-Poisson.h-0.50-coD1.00_1D-N000.txt` of set `A`. For more details, get the help of function `procWaveWTMM1d` by typing:

	> help procWaveWTMM1d
The function performing the main computation is stored in the file `myAnalyzeSeriesWaveWTMM`; type

	> help myAnalyzeSeriesWaveWTMM
for help and have a look on the code inside this file. Note that this method uses only the classical method to compute the Legendre transform (not the canonical one, like in _Lastwave_).

You can have a look in the test file, whose purpose is to use different implementations (but all using the main _WaveLab_ functions) for retrieving similar results.

*fraclab*

The directory `fraclab/` contains the main program for analyzing the signals with the WTMM method implemented in the [_FracLab_](FracLab) software.

First of all, you should be sure that the fraclab functions are known in the _Matlab_ environment. For this purpose, you have two solutions:

  1. either you install the library _FracLab_ in your environment; then, 
  you need to add the path of the directory _FracLab_ to the `PATH` variable 
  of _Matlab_, _e.g._ thanks to the command addpath (See the help of function 
  `procFracWTMM1d`, below).
  2. you use the 'stand-alone' files provided in the sub-directory 
  `fraclib/` of the directory `fraclab/`. 

In _Matlab_, move to the directory `fraclab/`. Then, in order to perform the estimation for all the series, you should type: 

	> procFracWTMM1d('all','all');
The processing of all the series takes about 25mn. However, you may want to try first over small dataset to begin, so type 
for example:

	> procFracWTMM1d('A',1); 
to perform the analysis only for the multifractal stored in file `Log-Poisson.h-0.50-coD1.00_1D-N000.txt` of set A. For more details, get the help of function `procFracWTMM1d` by typing:

	> help procFracWTMM1d
The function performing the main computation is stored in the file `myAnalyzeSeriesFracWTMM`; type

	> help myAnalyzeSeriesFracWTMM
for help and have a look on the code. Note that this method uses only the classical method to compute the Legendre transform.

We however recommand first to have a look on the `README` of the sub-directory `test/` and to perform the operations described therein.  The aim of the described operations is to help you to understand how _FracLab_ works, to see the results it provides and eventually see some of its limitations (not the canonical one, like in lastwave).

[LastWave]: http://www.cmap.polytechnique.fr/~bacry/LastWave/
[WaveLab850]: http://statweb.stanford.edu/~wavelab/
[FracLab]: https://project.inria.fr/fraclab/

**<a name="Reference"></a>References**

* Turiel A., Pérez-Vicente C.J., and Grazzini J. (2006): [**Numerical methods for the estimation of singularity spectra on sampled data: a comparative study**](http://www.sciencedirect.com/science/article/pii/S0021999105005565), _Journal of Computational Physics_, 216(1):362-390, doi: [10.1016/j.jcp.2005.12.004](http://dx.doi.org/10.1016/j.jcp.2005.12.004).
* Arneodo A., Audit B., Decoster N., Muzy J.-F. and Vaillant C. (2002): [**Wavelet based multifractal formalism: Application to DNA sequences, satellite images of the cloud structure and stock market data**](http://germain.its.maine.edu/~khalil/courses/MAT500/papers/arneodo_bookfractals_02.pdf), in Bunde A. et al. eds.: _The Science of Disasters: Climate Disruptions, Heart Attacks, and Market Crashes_, Springer Verlag, pp. 26-102. 
* Lévy Véhel J. and Tricot, C. (2004): [**On various multifractal spectra**](https://hal.inria.fr/inria-00559102/file/jlvctfinal.pdf), _Progress in Probability_, Fractal Geometry and Stochastics III, Birkhauser Verlag, 57, pp.23-42.
* Abry P., Gonçalves P., Lévy Véhel J. (2002): **Lois d’Echelle, Fractales et Ondelettes**, Hermès Sciences Publications.  
* Arneodo A., Bacry E., Jaffard S. and Muzy J.-F. (1998): [**Singularity spectrum of multifractal functions involving oscillating singularities**](http://link.springer.com/article/10.1007/BF02475987), _Journal of Fourier Analysis & Applications_, 4(2):159-174. doi:[10.1007/BF02475987](https:/doi.org/10.1007/BF02475987).
* Bacry E., Muzy J.F. & Arneodo A. (1993): [**Singularity spectrum of fractal signals from wavelet analysis: Exact results**](http://link.springer.com/article/10.1007/BF01053588), Journal of Statistical Physics, 70(3):635–674, doi:[10.1007/BF01053588](https:/doi.org/10.1007/BF01053588).
* Donoho D.L. (1993): [**Nonlinear wavelet methods for recovery of signals, images, and densities from noisy and incomplete data**](https://statistics.stanford.edu/sites/default/files/EFS%20NSF%20437.pdf), in Daubechies I. (ed.): _Different Perspectives on Wavelets_, American Mathematical Society, pp.173-205.
 
