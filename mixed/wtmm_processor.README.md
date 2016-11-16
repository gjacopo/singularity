Program: `wtmm_processor`

## Warning:

The wtmm-related library implements the Wavelet Transform Modulus Maxima analysis
of 1D signals only. The WTMM approach consists essentially in:
    1. locate the local maxima of the absolute value of the Wavelet Transform,
    2. track maxima lines for increasing wavelet scale,
    3. compute the multifractal exponents (singularity spectrum) trhough Legendre
     Transform.

The parts 1 and 2 can be discussed for the choice of some parameters like: the 
wavelet filter, the minimum and maximum scales of analysis. At the same time, the 
partition function is computed, so this implies the choice of appropriate moments.
The part 3 implies the choice of a method for estimating the Legendre Transform:
direct approach or canonical approach. Both are available in our code. The last one 
has been implemented in LastWave: it fully inspired the version we implemented in our
code. 

Note that the code is not optimized in term of time consuming (in particular, wavelet 
convolution doesn't use FFT: i wanted to reduce uncertainty due to supplementary 
numerical computation; it's probablyy stupid...) 

In the following, we provide a description of the stand-alone program `wtmm_processor` 
for running WTMM. However, the comments below also hold for the shared program 
`multifractal_processor`.

Anyway, have a look on the code itself: the functions are fully described.

## Required files:
* `utils.c` : useful functions for allocations, storage, registration... some of
    them are inspired (=copied) by the msm library.
* `parse_args_wtmm.c` : list of functions for parsing the parameters of the WTMM
    method.
* `extrema.c` : list of functions for computing the continuous wavelet transform 
    and the extrema computation through maxima lines representation
    (i.e. the so-called WTMM): this is the main contribution of the method
    and corresponds to parts 1 and 2 of the WTMM scheme.
* `approx_legendre.c` : list of functions for computing the multifractal exponents
    through the Legendre Transform: this is the other important part; two methods 
    are implemented: direct and canonical (note that in _LastWave_, only the canonical 
    one was implemented: this directly inspired the code written in this file); this 
    corresponds to the part 3 of the WTMM scheme; every time it's possible, the 
    references to the original _LastWave_ codes are available (see functions `canonMeanCompute` 
    and `canonSpecCompute`).
* `multifractal_wtmm.c` : some important functions for launching the computations. 
* `const_wtmm.h` : stores some pseudo-constant useful for default parametrization
* `variables_wtmm.h` : contains all the global variables used in the program.

## Compiling command :
While being in the directory multifractal/, type:

    make wtmm

See file `../multifractal_processor.README.md` for further details.

## Run-time parameters 

You can get a brief description of the parameters by typing `wtmm_processor -h`.
In the following you will see a more detailed description of this help.

The parametrization of the program is similar to that of multifractal_generator: 
the function for parsing arguments to the program is based on the functions 
parse_* of the msm-related library; see file parse_args_wtmm.c

In the following, we describe the input parameters of the WTMM method and the 
way they have to be passed to the program wtmm_processor. Note that, unless
we point out on it, THE WTMM PARAMETERS HAVE TO BE PASSED TO THE PROGRAM
multifractal_processor WITH THE SAME ARGUMENTS. 

    Usage: ./wtmm_processor [-f file_name] [-N #series] [-txt] [-canon] 
        [-wav wav_index] [-der deriv_order] [-sc0 scale_0] [-scmax scale_max] 
        [-dsc scale_step] [-tsc time_scale] [-rsc ratio] [-q min_moment max_moment] 
        [-dq moment_step] [-regsc sc_min sc_max] [-cdh shift_spectrum] [-tauq] [-ver]

### General input parameters 
* `-f` : The name of the input file to be processed. It's a generic (base) name if the
    number _NSERIES_ (option `-N`) is >1 or if it's a binary file. It will be completed 
    with the 'usual' strings `-N???`, and also `.dat` or `.txt` according to the type  
    of the file: see option `-txt` for flag _BINARY_. Otherwise, if it represents a 
    unique text file to be processed, the string won't be completed and its size is 
    automatically computed.
* `-txt` : Flag. If enabled, the input data are supposed to be in text format,
    otherwise they are in binary format. Default: _DISABLED_ 1

NOTE: THE FLAG `-f` AND `-txt` ARE NOT USED IN THE SHARED PROGRAM

* `-N` : Number of series to be processed. Default: 1

NOTE: THIS FLAG IS PROVIDED BY THE REGULAR FLAG `-N` OF `multifractal_generator`
IN THE SHARED PROGRAM

* `-ver` : Flag. If enabled, the program shows lots of (verbose) information,
specially for multifractal analysis. Default: _DISABLED_

NOTE: THIS FLAG IS PROVIDED BY THE REGULAR FLAG `-ver` OF `multifractal_generator`
IN THE SHARED PROGRAM

### WTMM parameters 
* `-canon` : Flag. If enabled, the method used in WTMM scheme to approximate the
    Legendre Transform is the canonical approach. Default:  _DISABLED_ (direct method)
    The direct estimation of multifractal spectrum is realized through the 
    estimation of the exponent tau(q) from the partition function. The canonical
    doesn't need these exponents (even if it provides it), it estimates the spectrum
    by runing a regression over some averaged exponents.
    [See the difference between the function `directSpecCompute` and the functions 
    `canonMeanCompute and `canonSpecCompute` in file `approx_legendre.c`. See also the 
    call to these functions in function estima_Dh_WTMM of file multifractal_wtmm.c.]
* `-wav` : Wavelet to be used in WTMM scheme, either:     
    + -1: Lorentzian,
    + 0: Morlet,
    + 1: Gaussian.
  
    Default:  1
    [See the function `define_wavelet` and `comp_coefficient` of `extrema.c` for description.]
* `-der` : Derivative order, only when the wavelet used in WTMM scheme is the gaussian.
  You can choose till the 5th derivative. Default:  2

NOTE: THE FLAGS `-wav` AND `-der` ARE REPLACED BY `-wavwtmm` AND `-derwtmm` IN THE SHARED 
PROGRAM RESP.

* `-sc0` : Initial scale of analysis in WTMM scheme. Default:  2
* `-scmax` : Maximal scale of analysis in WTMM scheme. By default, if it is not specified,
  it is computed automatically according to the lenght of the signal (see options 
  `-tsc` and `-scrat` below).
  [See the default definition of the maximal (when it is not provided) in function 
  `checkWTMMparameters` of file `multifractal_wtmm.c`.]
* `-dsc` : Scale step of analysis in WTMM scheme. Default:  1.0353. 
  The scales of analysis will range from sc0 to scale_max by a multiplicative step 
  scale_step.

NOTE: THE FLAGS `-sc0`, `-scmax` and `-dsc` ARE REPLACED BY `-sc0wtmm`, `-scmaxwtmm`
AND `-dscwtmm` IN THE SHARED PROGRAM RESP.

* `-tsc` : Factor of convolution range in WTMM scheme, i.e. defines the size of the wavelet 
  filter:

        [wavelet box] = #{-time_scale*nsc,...,time_scale*nsc}.
  Default:  5
* `-rsc` : Ratio between maximum scale and signal length. Default:  0.0625
  Naturarly, both parameters time_scale and ratio are taken into account for the 
  computation of the default maximal scale of analysis. 
  [See the function `checkWTMMparameters` of file `multifractal_wtmm.c`.]
* `-q` : Born of the range of moments used to perform the estimation of the multifractal
  exponents in WTMM scheme. Default:  [-5,5]
* `-dq` : Moment step of the WTMM analysis. If moment_step is choosed as 0, then a list
  of (irregularly spaced) moments will be used instead of the range of moments
  `[min_moment, max_moment]`.
  Default:  0.2
  [For an irregularly spaced list of moments, you must define it "in hard" in the  
  program. See variable `qDefArray` in file `variables_wtmm.h`]

NOTE: THE FLAGS `q` AND `-dq` ARE RESP. REPLACED BY `-qwtmm` and `-dqwtmm` IN THE SHARED 
PROGRAM

* `-regsc` : Born of the range of scales used to perform the estimation (through
  regression) of the multifractal exponents in WTMM scheme. Default:  [3,20]
  Only the scales in the range `[sc_min, sc_max]` will be used in the regression
  for estimating multifractal exponents.
* `-cdh` : Constant to add to the values of the WTMM spectrum when storing it.
  Default:  0
  Put `shift_spectrum=-1` to compare the results with spectra computed from msm-related 
  methods.
* `-tauq` : Flag. If enabled, the multifractal exponent computed by the WTMM scheme
  will be saved. Default:  _DISABLED_

NOTE: IN THE SHARED PROGRAM, THE WTMM METHOD FOR ESTIMATING MULTIFRACTAL EXPONENTS IS
APPLIED ON GENERATED SIGNALS ONLY IF THE FLAG `-analysis` IS ACTIVATED.

## Typical outputs

Remember that `$$$` the general basename associated to the input files (See option '-f'
below in the case of stand-alone program). 
When these inputs are produced by `multifractal_generator`, it consists of the string 
identifying the type of multifractal ("Log-Poisson", "Log-Normal" and "Log-Levi"), plus 
the tags and values for the associated free parameters, and a final tag specifying if 
the data are onedimensional (1D). 

* `Dh_wtmm_$$$` : singularity spectrum obtained from the WTMM method. The first column is the 
value of the singularity `h`, the second stands for the singularity spectrum `D(h)`. 

* `Tq_wtmm_$$$` : multifractal exponents `tau(q)` obtained from the WTMM method. They are 
defined from the partition function, and they are supposed to be related to the exponents 
`xi(p)` of the structure function by: `tau(q)=xi(q)-1`. They are the exact Legendre 
transform of the singularity spectrum `D(h)`. Thus, in the direct approach, they are used
for estimating the singularity spectrum `D(h)`. The first column is the moment value (taken 
in the list `qArray`), the second stands for the values of the exponents. This file is 
generated only if the flag `-tauq` is activated. 

## Typical execution examples

    ./wtmm_processor.exe -f Log-Poisson.h-0.50-coD1.00_1D -N 10 -der 2 -tsc 6. -regsc 4. 100. \
        -cdh -1. -tauq -ver 
performs WTMM over the files `Log-Poisson.h-0.50-coD1.00_1D-N00X.dat`, with `X={0,..,9}`,
with direct approach, using the 2nd derivative of gaussian; the scales ranging from 4
to 100 will be used to estimate the multifractal exponents; the outputs will consist 
in two files called `Dh_wtmm_Log-Poisson.h-0.50-coD1.00_1D` and 
`Tq_wtmm_Log-Poisson.h-0.50-coD1.00_1D`; except `time_scale=6`, the other parameters are 
the default ones. The spectrum is represented with a shift of -1.

    ./wtmm_processor.exe -f log-normal.txt -txt -der 3 -q -2. 8. -dq 0.5
performs the WTMM analysis of a text file called `log-normal.txt` with the 3rd derivative
of gaussian; the moments used for estimating the partition function will range from -2
to 8 with a step of 0.5; the output will be `Dh_wtmm_log-normal.txt`

    ./wtmm_processor.exe -f Log-Poisson.h-0.50-coD1.00_1D -txt -dq 0.
performs the WTMM analysis over the unique file `Log-Poisson.h-0.50-coD1.00_1D-N000.dat`;
in this case, as moment_step=0, the moments used for the estimation of the partition
function are those stored in the defaut list `qDefArray` (see `variables_wtmm.h`).

    ./multifractal_processor.exe -dermode 1 -N 1 -dim 16000 -type 1 -analysis -derwtmm 2 -tsc 6.\
        -regsc 4. 100. -cdh -1. -tauq -ver 
is similar to the first example above as regards WTMM analysis

    ./multifractal_processor.exe -N 10 -dim 16000 -analysis -derwtmm 3 -tsc 6. \
        -qwtmm -2. 8. -dqwtmm 0.5 -ver 
performs an analysis on Log-Poisson signals with similar parameters for the WTMM analysis
as in the second example above.

    ./multifractal_processor.exe -N 1 -dim 16000 -analysis -derwtmm 2 -tsc 6. \
        -dqwtmm 0 -regsc 4. 100. -ver 
can be used to perform the analysis of a default Log-Poisson signal with the list of 
moments:

    [-4.,-3.6,-3.2,-3.,-2.8,-2.6,-2.4,-2.2,-2.,-1.8,-1.6,-1.4,-1.2,-1.1,-1.,-0.9,-0.8,
    -0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,
    0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.,1.05,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,
    2.,2.3,2.6,3.,3.5,4.,5.,6.,7.,8.]
    
and with the range of scales `[4.,100.]` for regression run.

Version: Dekembriou 2th, 2004 - 24 degrees... 
