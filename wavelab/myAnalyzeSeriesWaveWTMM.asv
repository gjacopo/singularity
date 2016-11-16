% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% MYANALYZEWAVESERIESWTMM - WTMM analysis for multifractal spectrum estimation.
% Computes the WT and the extrema representation of a family of time series, 
% and then calculates the corresponding multifractal spectra.
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% The main processing parts of the code make an intensive use of the WaveLab 
% library of Matlab functions, available at:
%            http://www-stat.stanford.edu/~wavelab/
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Usage
%     [ h, d, tau, Z ] = AnalyzeSeriesWaveWTMM( storedir, savedir, base, ...
%                  binary, nseries, N, Q, R, wavelet, scale0, nvoice, noct )
% Input
%    storedir name of the directory where the data will be loaded
%    savedir  ____ __ ___ _________ _____ ___ results will be stored
%    base     generic name of the files storing the input (real-valued) signals  
%    binary   Flag 0/1 for input format: 0 for text, 1 for binary 
%    nseries  number of series to analyze
%    N        size of the time series; if unknown, put N<=0 
%    Q        list of exponents
%    R        vector 2x1 of min and max scales to fit line by
%    wavelet  analyzing wavelet among 'Gauss','DerGauss','Sombrero'
%    scale0   minimum scale of analysis (def=2)
%    nvoice   number of voices/octave (def=10)
%    noct     no of octaves discarded by the analysis (def=2)
% Output
%    h        vector 1 x nscale of h: Singularity Exponents
%    d        vector 1 x nscale of d(h): Legendre Spectrum
%    tau      vector 1 x nq of tau(q): Multifractal Exponents
%    Z        matrix nexp x nscale of z(q,a): Partition Function
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [ h, dh, tau, Z ] = myAnalyzeSeriesWaveWTMM( storedir, savedir, base, ...
        binary, nseries, N, Q, R, wavelet, scale0, nvoice, noct )
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Internal RELEVANT parameters
REVERSE=1;
par=100000000;
MYIMPLEMENT=1;
SKEL=0; % flag for the extraction of the ridges
        % Fixing the variable SKEL to 0 may increase dramatically the
        % computation time

HMIN=-2.; HMAX=1.5;

% Internal parameters/flags: you can eventually pass them as parameters to
% the function, but they are not relevant...
flagDISPLAY=1; % display the results or not
precision='double'; % format of the input data

if(exist('R') ~= 1), R=0; 
elseif (size(R,1)== 1 & size(R,2)== 1), R=0; end;

% check parameters
if(exist('scale0') ~= 1), scale0=2; end;
if(exist('nvoice') ~= 1), nvoice=10; end;
if(exist('noct') ~= 1), noct=2; end;

cut0=scale0;;  

% Eventually determine the length of the time series if N hasnt been set 
if (N <= 0)
    N = dimfileseries( base );
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
     if REVERSE, % miror the signal
       x(N+1:2*N) = reverse(x);
        noctave=noctave+1;
    end
    nscale  = nvoice * noctave;
    
    % Computes the Continuous Wavelet Transform of the signal x
    % wt = RWT( x, nvoice, wavelet, nscale, scale0);	
    wt = CWT( x, nvoice, wavelet, noct, scale0 );
    % Inputs
    %    x        signal, dyadic length n=2^J, real-valued
    %    nvoice   number of voices/octave
    %    wavelet  string 'Gauss', 'DerGauss','Sombrero', 'Morlet'
    %    nscale   Default=2
    %    scale    Default=4
    % Outputs
    %    wt      matrix n by nscale where nscale = nvoice .* noctave
 
     if REVERSE
         % Resize the matrix of WT coefficients
         [n,nscale] = size(wt);
         wt = wt(:, cut0*nvoice + 1 : nscale - nvoice);
         [n,nscale] = size(wt);
         scale = 2.^(log2(scale0) + (noctave-cut0-1)*(1:nscale)/nscale);
     else
         scale = 2.^(log2(scale0) + noctave*(1:nscale)/nscale);
      end;
      
    if (flagDISPLAY && nseries==1)
        % display CWT   
        scaling='Individual'; color='hot'; % opt='lin';
        figure; ImageCWT( wt, scaling, color );
        % IMAGECWT -- Image of Continuous Wavelet Transform
        title('CWT')
    end;
    
    % Build Maxima Map relating the positions of the WTMM
    maxmap = MM_RWT( wt, par );
    % Returns maxmap: binary array indicating presence of max or not
    % Note: ne pas utiliser maxmap = WTMM( wt );

    if (flagDISPLAY && nseries==1)
        % display maxmap, i.e. the maxima lines
        figure; ImageWTMM( fliplr(maxmap) )
        % IMAGEWTMM -- Maxima of Continuous Wavelet Transform
    end;
      
    % Build Thermodynamic Partition Function
    % Computes the Fractal Partition Function based on wavelet modulus maxima
    %             z(q,a) = sum_i( |CWT(a,b(i))|^q )
    % where b = (b(i)) is a list of wavelet transform maxima.  
    if SKEL
        % probably better, but also longer
        part = myCalcThermoPartition( wt, maxmap, Q );
    else
        part = CalcThermoPartition( wt, maxmap, Q );
        % Note: This function is similar to the function FRACPARTITION of 
        % WaveLab802/Fractal/:
        %    part = FracPartition( wt, maxmap, Q );       
    end;   
   % Returns part: matrix nexp by nscale of z(q,a), where nexp=length(Q) and
   % nscale=size(wt,2)

    % Update the statistics of the overall series: simply add the partition
    % functions:
    if (ns==1) Z = part;
    else   Z = Z + part;
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
end    
if isempty(R) | R==0, 
    loscale = scale(1); 
    hiscale = scale(nscale);
end;

% Computes the multifractal exponents: Calculate Moment Generating Function
%    tau(q) = Slope [ log(z(q,a))  versus log(a) ]
if MYIMPLEMENT
    tau = myCalcMomentGenFun( Z, scale, loscale, hiscale );
    % note: tau=-tau is done inside this function
else
    tau = FracScalExp( Z, scale, loscale, hiscale );
    % Note: This function is identical to the function FRACSCALEXP of 
    % WaveLab802/Continuous/:
    %         tau = CalcMomentGenFun( Z, scale, loscale, hiscale );
    % Note that we have to invert the exponents:
    tau = -tau;
end
% Returns tau: vector 1 by nq, where nq=size(z,1) is the no of moments
    
% Calculate Spectrum of Local Scaling Exponents
% Compute the multifractal spectrum:
%        h = \frac{d\tau}{dq}
%        D(h) = q h -tau(q) 
if MYIMPLEMENT
  [ h, dh ] = myCalcSingSpectrum ( tau, Q );
else
    h = linspace( HMIN, HMAX, 101 );
    dh = FracSingSpect( tau, Q', h );
    % Note: This function is similar to the function CALCFRACSPECTRUM of 
    % WaveLab802/Continuous/
    %    dh = CalcFracSpectrum(tau,Q',alpha);
    % ... modulo Q/2!
    % Get the dimensions >0
    i = find(dh >= 0);
    dh = dh(i);    h = h(i);
end
% Returns h: the list of singularity exponents values and dh: the singularity
% spectrum

if (flagDISPLAY)
  figure; plot( Q, tau );
  text = sprintf('%s - Multifractal Exponents',base); 
  title(text);
  xlabel('q');	ylabel('Tau(q)');
  % PlotMomentGenFun( tau, Q );
 
  figure; plot( h, dh );
  text = sprintf('%s - Multifractal Spectrum',base);
  title(text);
  xlabel('h');	ylabel('D(h)');
  % PlotMultiSpectrum( h, dh );
end

savefilevar( savedir, 'Tq_wave_', base, Q, tau );
savefilevar( savedir, 'Dh_wave_', base, h, dh );


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

