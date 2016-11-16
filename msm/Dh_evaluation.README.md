Program: Dh_evaluation

## Compiling command (using cc):

	gcc -lm -I<path> -o Dh_evaluation.exe Dh_evaluation.c

where the chain "<path>" must be substituted by the path to reach
the directory with the libraries required by this program.

## Required libraries

	struct_def.c
	operaciones.c
	tensor.c
	FFT.c
	derivacion.c
	multifractal.c

## Run-time parameters 

You can get a brief description of the parameters by typing:

	Dh_evaluation -h
In the following you will see a more detailed description of this help.

There are three types of run-time variables: boolean, integer and floating 
point. 

Boolean variables can handle two states only, _ENABLED_ or _DISABLED_. By
default, they are all _DISABLED_; they switch to the state _ENABLED_ when the
corresponding run-time parameter is specified; this parameter does not allow
accompanying numerical factors.

Integer and floating point variables must be accompanied by a numerical value
of the appropriate type, which specifies its value. Those variables have
default values, which can be known when the "-h" run-time parameter is passed.

Besides, some variables are sons of boolean variables, that is, if a value
is specified for one of these variables, automatically the state of the mother
boolean variable becomes _ENABLED_, even if it was not already specified. This
is a convenient behaviour for some cases, in which the activation of a given
processing mode opens the possibility of modifying some parameters specific to
that mode; with this automatic enabling feature you can directly modify the
parameters specific for that mode and you do not need to explicitly enable it
(what sometimes you can forget): once a son parameter is modified the mode is
enabled. 

	Usage: Dh_evaluation [-ver] [-fromDh] [-geomap] [-typemap] [-N #series] 
		[-dim length] [-d_space dimension] [-type mult_type] [-hinf min_sing] [-Codinf
		min_sing_cod] [-h1 max_sing] [-mean sing_av] [-sigma disp_sing] [-alpha
		exponent] [-memory] [-wav wav_index] [-der deriv_order] [-hold] [-s0 scale_0]
		[-range scale_rng] [-Nbin #bins] [-Method ] [-LastWave] 

### General parameters

 * `-ver` : Flag. If enabled, the program shows lots of (verbose) information,
	specially for multifractal analysis. Default: _DISABLED_
 * `-fromDh` : Flag. If enabled, the program takes previously computed `D(h)` files
	and estimates the error from them.
	Default: _DISABLED_
 * `-geomap` : Flag. If enabled, the program tries to generate quality maps
	for each method, changing geometry but keeping the given MF type.
 	Default: _DISABLED_
 * `-typemap` : Flag. If enabled, the program tries to generate quality maps
	for each method, for fixed geometry (1x16384) and changing parameters in
	the given MF type.
 	Default: _DISABLED_
 * `-N` : Number of series to be processed. Default: 1
 * `-dim` : Size of series to be processed. Default: 512
 * `-d_space` : Dimension of the embedding space. Default: 1
 * `-type` : Type of multifractal to be generated.
   	0: Log-Poisson
   	1: Log-Normal
   	2: Log-Levi
   	3: Binomial
 	Default: 0
 * `-hinf` : Most singular exponent. Valid for log-Poisson and binomials. Default: -0.50
 * `-Codinf` : Most singular codimension. Valid for log-Poisson and binomials. Default: 1.00
 * `-h1` : Maximum singularity in binomial MFs. Default: 0.50
 * `-mean` : Singularity mean. Valid for log-Normal and log-Levi. Default: 0.50
 * `-sigma` : Singularity dispersion. Valid for log-Normal and log-Levi.
	Default: 1.00
 * `-alpha` : Valid for log-Levi only: log-Levi exponent. Default: 1.50

### Memory parameters
* `-memory` : Flag. If enabled, the program makes a better use of memory (at the
cost of longer processing times). Default: _DISABLED_

### Wavelet variables
* `-wav` : Wavelet to be used.
   	0: Gaussian
   	1: Lorentzian at exponent 0.5
   	2: Lorentzian
   	3: Lorentzian at 1.5.
	Default:  0 
* `-der` : Derivative order for the wavelet. Default:  0
* `-hold` : Flag. It enables function processing (in opposition to
	measure analysis). Default: _DISABLED_
* `-s0` : Initial scale for wavelets. Default:  1.000000
* `-Nbin` : Number of histogram bins. Default: 128
* `-Method` : Composite variable expressing the method or methods to be used. 
	Its value is `1*GH+2*GMWP+4*MOM+8*WTMM`. Default: 15

### WTMM parameters 
* `-LastWave` : Flag. If enabled, WTMM calculations are done by external calls to
 	LastWave. Default: _DISABLED_

## Keyboard inputs

First (and only) input: Path to the directory in which series are stored.

### File inputs

The program assumes that a collection of files with appropriate names will be
placed at the given directory. Names and formats must correspond to those
produced by the program `multifractal_generator``.

### Typical outputs

The four methods to evaluate the singularity spectrum can be invoked; each one
has an associated alphanumerical chain (namely, `"Mom", "WTMM", "GMWP" and
"GH"`) which will be attached to the corresponding output filename. Recall that
you can activate any number of these calculations methods by using the
run-time parameter `-Method`. In the present version, the program incorporates a
built-in version of WTMM. You can change this however to a call to an external
instalation of _LastWave_; in such a case, you should specify `-LastWave` as a
run-time flag. You may need to modify `LWCALL` in line 84 to accomodate to the 
name of the LastWave command in your system.

In any working mode you will be provided the singularity spectrum files
associated to the method used and the ensemble analyzed. The name of a
spectrum file is of the form _Dh_<Method>_<Ensemble>_ where <Method>
designates the validation method and <Ensemble> defines all the
characteristics of the analyzed ensemble. Spectrum files are four-column,
ascii-formatted files. The columns mean, by order, `h`, `D(h)`, errorbar and
theoretical `D(h)`. For methods which do not allow to evaluate errobars, the
third column is conventionally fixed to 0.2.

We will represent by "$$$" the general basename associated to the output
files. It consists of the name identifying the type of multifractal
("Log-Poisson", "Log-Normal" and "Log-Levi", in the present version), plus the
tags and values for the associated free parameters, and a final tag specifying
if the data are onedimensional (1D) or bidimensional (2D).

Modes `-geomap` and `-typemap` have an special effect. When you activate -geomap
you are suppose to evaluate a benchmark of ensembles, according to the
geometry defined by the Size Benchmark in the paper. The flag -typemap
evaluates spectra according to the prescription of Spectrum Benchmark of the
log-Poisson or the log-Normal types (no other are given so far). So, within
these two modes, in addition to the spectra for each ensemble in the benchmark
you will obtain files `GeoMap-$$$` and `TypeMap-$$$` with the estimation errors
for each element in the benchmark. The files are formatted in ascii.

Mode `-fromDh` has been introduced to ease calculations of GeoMaps and
TypeMaps. If you have already calculated the spectra and the files are in the
same directory, the flag `-fromDh` will process these files and obtain the maps
of errors from them. This allows fast re-calculations, and may be very
convenient to avoid re-starting the process if the program hungs.
