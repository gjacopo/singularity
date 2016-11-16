% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% MYANALYZESERIESFRACWTMM - WTMM analysis for Legendre Multifractal Spectrum 
% determination based on the FRACLAB library.
% Computes the WT and the extrema representation of a family of time series, 
% and then computes the corresponding multifractal spectra.
% See the functions MYCWTPART and MYCWTSPEC, MYCONTWTPART and MYCONTWTSPEC.
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% The main processing parts of the code make an intensive use of the Fraclab 
% library of Matlab functions, available at:
%            http://fractales.inria.fr/index.php?page=fraclab
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Usage
%     [ h, d, tau, Z ] = myAnalyzeSeriesFracWTMM( storedir, savedir, base, ...
%           binary, nseries, N, Q, R, scale0, nvoice, noct, mode )
% Input
%    storedir name of the directory where the data will be loaded
%    savedir  ____ __ ___ _________ _____ ___ results will be stored
%    base     generic name of the files storing the input (real-valued) signals  
%    binary   Flag 0/1 for input format: 0 for text, 1 for binary 
%    nseries  number of series to analyze
%    N        size of the time series; if unknown, put N<=0 
%    Q        list of exponents
%    R        vector 2x1 of min and max scales to fit line by
%    scale0   minimum scale of analysis (def=2)
%    nvoice   number of voices/octave (def=10)
%    noct     no of octaves discarded by the analysis (def=2)
%    mode     'L1' or 'L2'. Have a look on files MYCWTSPEC and MYCONTWTSPEC.
% Notes: 
% i) the total number of octaves used in the CWT is floor(log2(N))-noct.
% ii) here, we use only the mexican hat as the analyzing wavelet, and we 
% don't parameterize the minimal scale, the number and voices or octaves:
% it means still more parameters to adjust...
% Output
%    h        vector 1 x nscale of h: Singularity Exponents
%    d        vector 1 x nscale of d(h): Legendre Spectrum
%    tau      vector 1 x nq of tau(q): Multifractal Exponents
%    Z        matrix nexp x nscale of z(q,a): Partition Function
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [ h, dh, tau, Z ] = myAnalyzeSeriesFracWTMM( storedir, savedir, base,  ...
        binary, nseries, N, Q, R, scale0, nvoice, noct, mode )
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Internal parameters
% you can eventually pass these ones as parameters to the function, 
REVERSE=1;
% irelevant parameters
flagDISPLAY=1; % display the results or not
precision='double'; % format of the input data

if(exist('R') ~= 1), R=0; 
elseif (size(R,1)== 1 & size(R,2)== 1), R=0; end;

% check parameters
if(exist('scale0') ~= 1), scale0=2; end;
if(exist('nvoice') ~= 1), nvoice=10; end;
if(exist('noct') ~= 1), noct=2; end;
if(exist('mode') ~= 1), mode='L2'; 
elseif (~strcmp(mode,'L1') & ~strcmp(mode,'L2')), mode='L2'; end;

cut0=scale0;

% Eventually determine the length of the time series if N hasnt been set 
if (N <= 0)
    N = dimfileseries( storedir, base );
end;
% In the analysis scheme, only N data will be loaded from the different
% files, even if they contain more data.

if(flagDISPLAY), fprintf('\nAnalyzing the processes %s - set %s ...',base,storedir);
end;

for ns = 1:nseries
 
    %% Load the file
    % The signal x must be real-valued with dyadic length N=2^j
    [x, bb] = loadfileseries( storedir, base, ns-1, N, binary, precision );
    if(flagDISPLAY), fprintf('\n     processing file %s', bb); end;
         
    if (flagDISPLAY && nseries==1)
        figure; plot(x);
        text = sprintf('Signal %s',bb);
    	title(text);
    end

    noctave = floor(log2(N))-noct;
     if REVERSE
       x(N+1:2*N) = reverse(x);
        noctave=noctave+1;
    end
    nscale  = nvoice * noctave;
 
    %% Computes a continuous wavelet transform of a 1-D signal (real or complex
    % The scale operator is unitary with respect to the L1 norm (function #cwt#)
    % or L2 norm (function #contwt#). 
    % Two closed form wavelets are available in the original code: the Mexican 
    % Hat or the Morlet Wavelet (real or analytic). In the following, we will use 
    % the Mexican Hat only.
    %
    % Range of frequencies to perform the analysis. 
    % Eventually change these lines.
    fmin = 2^(-noctave);
    fmax = min( 0.5, 2^(-log2(scale0)));
    %
    % Note: if you choose one of the functions above instead of cwt1d, also
    % replace the functions mycwtpart and mycwtspec below by mycontwtpart and
    % mycontwtspec respectively.
     if (strcmp(mode,'L2')) 
         [ wt, scale ] = contwt( x, fmin, fmax, nscale, 0 );
         %  [ wt, scale ] = contwtmir(x,fmin,fmax,N,wvlt_length)
     elseif (strcmp(mode,'L1')) 
         [ wt, scale ] = cwt1d( x, fmin, fmax, nscale, 0 );
     end
     % Input parameters
    %   - x : Real or complex vector [1,nt] or [nt,1], signal to be analyzed.
    %   - fmin : real scalar in  [0,0.5], Lower frequency bound of the analysis. 
    %   - fmax :  real scalar [0,0.5], Upper frequency bound of the analysis.
    %   - nscale : positive integer, number of analyzing voices.  
    %   - wvlt_length  :       0: Mexican hat wavelet (real)
    % Output parameters
    %   - wt : Real or complex matrix [N,nt], coefficients of the wavelet transform.
    %   - scale : real vector [1,N], analyzed scales

     if REVERSE
         % Resize the matrix of WT coefficients
         wt = wt(cut0*nvoice+1 : nscale - nvoice,:);
         scale = scale(cut0*nvoice+1 : nscale - nvoice);
         [n,nscale] = size(wt);
     end
     
    %% Build Maxima Map and computes the Fractal Partition Function based on 
    %% wavelet modulus maxima
    %             z(q,a) = sum_i( |CWT(a,b(i))|^q )
    % where b = (b(i)) is a list of wavelet transform maxima.  
    if (strcmp(mode,'L2')) 
        part = mycontwtpart( wt, scale, Q, 1);
    elseif (strcmp(mode,'L1')) 
        part = mycwtpart( wt, scale, Q, 1);
    end
   % Input parameters   
    %   - wt : Wavelet coefficients of the CWT 
    %   - scale : Analyzed scale vector
    %   - Q : Exponents of the partition function
    %   - FindMax : FindMax=1 means tha we estimate the Legendre spectrum from 
    %     the local Maxima, and not from all coefficients.
    % Output parameters
    %   - part : Partition function

    % Update the statistics of the overall series: simply add the partition
    % functions:
    if (ns==1) Z = part;
    else      Z = Z + part;
    end;
    
end; % end of loop "for n=1:nseries"
 
%% Estimates the multifractal Legendre spectrum through regression on the
% partition function
% Determines the range of the regression
if R~=0,
    loscale = R(1);
    hiscale = R(2);
    window = (loscale <= scale) & (scale <= hiscale);	
    R     = find(window);
    if isempty(R), R=0; end;
end

%% Effective computation    
if (strcmp(mode,'L2')) 
    [ h, dh, tau ] = mycontwtspec( Z, scale, Q, R );
elseif (strcmp(mode,'L1')) 
    [ h, dh, tau ] = mycwtspec( Z, scale, Q, R );
end;
% Input parameters
%   - part : Partition function
%   - scale : Analyzed scale vector
%   - Q : Exponents of the partition function, i.e. list of moments
%   - ChooseReg : 0/1 flag or integer vector [1,N_reg], N_reg<=N_scale
%        ChooseReg = 0 : full scale range regression
%        ChooseReg = [n1 ... nN_reg] : scale indices for the linear
%        regression of the partition function.
% Output parameters
%   - h : Singularity support of the multifractal Legendre spectrum
%   - dh : Multifractal Legendre spectrum
%   - tau : Multifractal Exponents

if (flagDISPLAY)
  figure; plot( Q, tau );
  text = sprintf('%s - Multifractal Exponents',base);
	title(text);
	xlabel('q');	ylabel('Tau(q)');
  figure; plot( h, dh );
    text = sprintf('%s - Multifractal Spectrum',base);
	title(text);
	xlabel('h');	ylabel('D(h)');
end

savefilevar( savedir, 'Tq_frac_', base, Q, tau );
savefilevar( savedir, 'Dh_frac_', base, h, dh );



% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% LOADFILESERIES -- Loads binary file storing time series data
% Usage
%      [x,bb] = loadfileseries(storedir,base,index,N,binary,precision)
% Input
%    storedir    name of the directory where the data are stored
%    base        generic (base) name; it will be completed with the 'usual' 
%                strings '-N???', and also '.dat', '.fdat' or '.txt'
%                according to the type of file (see input binary) and the 
%                parameter format (see input precision)
%    index       index of the file to load; for instance, if index=1, then
%                the file <base>-N001.dat (ouo .fdat) will be load
%    N           given size of the input file (length of the time series) 
%    binary      format of the input file: 0 (text) or 1  (binary)
%    precision   storage of the times series: 'float' on 32 bits or 'double'
%                on 64 bits; default: 'double'
% Output
%    x           input signal loaded
%    bb          full name of the file (without the directory)
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [x,bb] = loadfileseries( storedir, base, index, N, binary, precision )
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if (nargin < 4 || exist('precision') ~= 1 || exist('binary') ~= 1)
    binary=1;
    precision='double';
elseif (nargin < 5 || exist('precision') ~= 1)
    precision='double';
end

if binary
    if (strcmp(precision,'float'))
        bb = sprintf('%s-N%03d.fdat',base,index);
    else
        % all other cases 
        bb = sprintf('%s-N%03d.dat',base,index);
        % for eventual bad strings
        precision='double';
    end        
    nombre = sprintf('%s/%s', storedir,bb);
    fid = fopen( nombre, 'rb');
else
    bb = sprintf('%s-N%03d.txt',base,index);
    nombre = sprintf('%s/%s', storedir,bb);
    fid = fopen( nombre, 'r');
end

% read the data and store in x
if binary
    x = fread( fid, N, precision );
else    
    if (strcmp(precision,'float'))
        x = fscanf(fid, '%f', N);
    else
        x = fscanf(fid, '%lf', N);
    end
end

fclose(fid);


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% DIMFILESERIES -- Determines the size of a file by reading the associated
% header file
% Usage
%      N = loadfileseries(storedir, base)
% Input
%    See function LOADFILESERIES
% Output
%    N     size of the input file (length of the time series) 
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function N = dimfileseries( base )
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Define the file to open
nombre = sprintf( '%s-N000', base );
ll = length( nombre );
nombreh = sprintf( '%s.hdr', nombre );

% Open the header file
fid = fopen( nombreh, 'r');

% Read the first line, which should correspond to the name of the file
nome = fgetl( fid );

if (~ischar(nome) || ~strncmp(nome,nombre,ll))    
    fprintf(['File name:\n\t %s \nin header file\n\t %s',... 
             '\ndoesn''t correspond to a known file'],...
             nome, nombreh );    
     N=-1;    
else    
    % Read the size of the time series
    N = fgetl( fid );
    % Convert the string to a number 
    N = str2num( N );
end

fclose(fid);

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% SAVEFILEVAR -- Store multifractal variables (exponents and/or spectrum) in
% a text file
% Usage
%      savefilevar(savedir, pref, base, x, y)
% Input
%    savedir     name of the directory where to save the results
%    pref        prefix to add to the name of the output
%    base        generic (base) name
% The data will be saved in 'double' format.
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function savefilevar( savedir, pref, base, x, y )
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

nombre = sprintf('%s/%s%s',savedir,pref,base);
z = [x(:), y(:)];

fid = fopen( nombre, 'w');

% save the data 
fprintf( fid, '%f %f\n', z' );
 
fclose(fid);

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function r = reverse(x) 
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  r = x(length(x):-1:1);

