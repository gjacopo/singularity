Program: `multifractal_generator`

## Compiling command (using cc):

    gcc -lm -I<path> -o multifractal_generator.exe multifractal_generator.c

where the chain _<path>_ must be substituted by the path to reach
the directory with the libraries required by this program.

## Required libraries

    struct_def.c
    operaciones.c
    tensor.c
    FFT.c
    derivacion.c
    multifractal.c
    graficos.c

## Run-time parameters 

You can get a brief description of the parameters by typing
`multifractal_generator.exe -h`. In the following you will see a more detailed
description of this help.

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

    Usage: multifractal_generator [-ver] [-N #series] [-dim size] 
        [-outres #res_levels] [-invtrans] [-d_space dimension] [-type mult_type] 
        [-hinf min_sing] [-Codinf min_sing_cod] [-h1 max_sing] [-mean sing_av] 
        [-sigma disp_sing] [-max_disp #sigmas] [-alpha exponent] [-memory] 
        [-dermode mode] [-wav wav_index] [-der deriv_order] [-hold] [-s0 scale_0] 
        [-wv_basis choice] [-der_wv_basis order] [-Nbin #bins]

### General parameters

* `-ver` : Boolean variable. When enabled, the program dumps more verbose
    information on standard output. Also, it produces printable representations
    of each generated series. WARNING: you can fast overcrow your hard disk.
* `-N `: Integer variable. Number of series to be generated. Default: 1
* `-dim` : Integer variable. Length for 1D series and linear size for 2D
    signals. The program approximates this value to the lowest integer greater
    than this and power of 2, for dyadic division purposes. Default: 512.
* `-outres` : Integer variable. Number of resolution levels to be ignored (by
    block reduction) in order to smoothen the generated series. Default: 0
* `-invtrans` : Boolean variable. If enabled, it forces the output data to be
    generated in a translationally invariant way. This implies the cancellation
    of the exponent `\tau_1`. Log-Poisson multifractals are constructed in a
    scale-invariant fashion, but this is not true for log-Normal and log-Levy;
    for these two types scale invariance is achieved by modifying the mean
    singularity value accordingly. In addition, if `-analysis` is specified the
    program does not try to compute the systematic shift produced by the lack of
    scale invariance in the derivative or the wavelets-on-measures methods. 
* ` -d_space` : Integer variable. Dimension of the generated signals. Default: 1
* `-type` : Integer variable. Type of multifractal to be generated.
        + 0: Log-Poisson
        + 1: Log-Normal
        + 2: Log-Levi
        + 3: Binomial
    Default:0
* `-hinf` : Float variable. Minimum singularity in log-Poisson MFs and binomial MFs. 
    Default: -0.5
* `-Codinf` : Float variable. Most singular codimension. Valid for log-Poissons and binomials. 
    Default: 1.00
* `-h1` : Maximum singularity in binomial MFs. Default: 0.50
* `-mean` : Float variable. Value for the singularity mean (_i.e._, maximum
    dimensionality). This parameter only affects the generation of log-Normal and
    log-Levi signals. Default: 0.50
* `-sigma` : Float variable. Singularity dispersion. It is valid for log-Normal
    and log-Levi multifractals only; in log-Normal it is equivalent to the
    concept of standard deviation but for final singularities; in log-Levi it is
    the closest equivalent. Default: 1.00
* `-max_disp` : Float variable. Singularity range, expressed in sigmas; it gives
    the effective range of explored singularities, which needs to be
    truncated: singularities can be observed from the mean minus this value times
    the sigma to the mean plus this value times the sigma. Take into account that
    even if the quasi-theoretical singularity spectrum seems correct you may be
    exhausting the numerical capabilities of the computer if this range is too
    wide (see below). This parameter only affects MFs of the types log-Normal and 
    log-Levi. Default: 5.00
* `-alpha` : Float variable, valid only if log-Levi is specified: log-Levi
    exponent, which controls the shape of the singularity exponent. For an
    exponent of 2 log-Normal distributions are retrived. Its extrem values are 0
    and 2. Default: 1.50

### Memory parameters
* `-memory` : Boolean variable. If enabled, the program uses FFFT1D instead of
    FFT1D when Fourier transforms are calculated (in convolutions, filters,
    etc). The difference between this two routines is not small: FFT requires
    the length of the series to be a power of 2, and if this is not the case,
    the program will be zero-padding the series up to the closest larger power
    of 2. On the other hand, FFFT does not require the dimension to be a power
    of 2, but the speed of the routine is dependent on the size of the prime
    factors in the decomposition of lenght. For lenghts which are powers of 2 or
    3, the speed of this routine is only slightly slower than that of FFT, and
    you save memory as the program does not perform zero-padding. This variable
    is not of great use in the present version of the program, as the sizes are
    always powers of 2 for other reasons; FFT is used by default.

### Derivative variables
* `-dermode` : Integer variable. Method to compute derivatives
        * 0: Half-pixel derivative kernel. It induces errors and non-local
        effects. However this mode is chosen by default for compatibility with
        other programs using this library; this fact must be taken into account
        when running this program.
        * 1: Direct one-pixel increment. This kernel coincides with the intuitive
        definition of derivative on series, and does not induce errors nor
        non-local effects. This sould be the mode of choice for the present
        program.
        * 2: Interpolation of derivatives. Noisy and experimental. Please ignore.

    Default:0

### Wavelet variables (used for compute singularity exponents)
* `-wav` : Integer variable. Wavelet to be used. 
        + 0: Gaussian 
        + 1: Lorentzian at exponent 0.5
        + 2: Lorentzian
        + 3: Lorentzian at exponent 1.5.
    (Recall that the slow the wavelet decays at infinity, the narrowest is the
    range of detected singularities but the best you can detect the MSM).
    Default:  0
* `-der` : Integer variable. Derivative order for the wavelet. Default: 0. If you
    choose to process signals instead of measures you should, at least,
    introduce on order of derivative.
* `-hold` : Boolean variable. When enabled functions are processed in a Holder
    scheme (in opposition to measure analysis)
* `-s0` : Float variable. Initial scale for wavelets. Default:  1.000000

### Parameters defining the representation basis
* `-wv_basis` : Integer variable. Wavelet of choice for the basis.
        + 0: Gaussian wavelet
        + 1: Lorentzian wavelet
        + 2: Diagonal Haar
    Default: 0
* `-der_wv_basis` : Integer variable. It indicates how many derivatives must be
    applied to the basiswavelet. Default: 2. Recall that this value must be
    greater than zero if you want to have a real wavelet basis (although this
    condition only does not guarantee it is a basis).
* `-Nbin` : Integer variable. Number of histogram bins. Default: 128

## Inputs

None

## Outputs

We will represent by `$$$` the general basename associated to the output
files. It consists of a root name identifying the type of multifractal
("Log-Poisson", "Log-Normal", "Log-Levi", and "Binomial"), tags and values for
the associated free parameters, and a final tag specifying if the data are
onedimensional (1D) or bidimensional (2D).

* `$$$-N?????.hdr` : Header file. The `?????` refers to a number between 00000 and
99999 which identifies that particular realization in a sequence up to 100000
(you fix the number of output series with `-N`). The file contains the name of
the associated data file and the linear length L of its records; that is, the
data file will contain L records if it is 1D and LxL if it is 2D.

* `$$$-N?????.dat` : Data file. It contains the data in raw binary double
format. You can substitute the routines `graba` and `lee` in the program with
the included "float" counterparts to produce and read float data. Take into
account that if you enable the `-analysis` run-time parameter the program will
need to re-read the data it has already generated and for that reason both
`graba` and `lee` routines have to be modified in the same way.

* `$$$-N?????.txt, $$$-N?????.gif` : If you have verbosed the results (`-ver`), the
program will generate a _plotable.txt_ file (for 1D series) or a GIF image (2D data). 

## Typical execution examples

    multifractal_generator -ver -dermode 1 -analysis -N 10 -dim 4000

generates 10 log-Poisson (with default choice of parameters) 1D series of 4096 elements and 
their associated plottable .txt files, and produces the different singularity 
spectrum estimates; it forces the program to derivate using 1-pixel increments

    multifractal_generator -ver -d_space 2 -analysis -type 1 -mean -0.25

generates one 512x512 2D log-Normal signal (`sigma=1, mean=-0.25`) and a
graphical representation, and its analysis. The choice of derivatives with 
1-pixel increments is not required as this is less influencing on 2D data

    multifractal_generator -dermode 1 -analysis -type 2 -invtrans -N 10 -dim 10000

generates 10 log-Levi 1D series of 16384 elements, with `alpha=1.5`, `sigma=1` and 
the mean adjusted by the requirement of translational invariance (`mean=0.33`), and 
analizes them.
