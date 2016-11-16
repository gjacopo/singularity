% TESTWAVEFUNCTIONS - Performs the multifractal analysis of some determined
% test signals with the original WaveLab code and some other derived code.
% Usage
%       testWaveFunctions(choice)
% Input parameter:
%  - choice : integer flag in {0,1,...8}
%        choice=1 : Log-Poisson.h-0.50-coD1.00_1D-N000.txt of set A
%        choice=2 : Log-Normal.mean0.50-sigma1.00_1D-N000.txt of set A
%        choice=3 : Log-Poisson.h-0.50-coD1.00_1D-N000.txt of set B
%        choice=4 : Log-Poisson.h-0.25-coD0.75_1D-N000.txt of set B
%        choice=5 : Log-Normal.mean0.50-sigma1.00_1D-N000.txt of set B
%        choice=6 : Log-Normal.mean0.25-sigma0.33_1D-N000.txt of set B
%        choice=7 : Log-Poisson.h-0.50-coD1.00_1D-N000.txt of set C
%        choice=8 : Log-Normal.mean0.50-sigma1.00_1D-N000.txt of set C
%        choice=0 : devilStaircase.txt
function testWaveFunctions(choice)

Nequiv={'devilStaircase.txt' ...
    'Log-Poisson.h-0.50-coD1.00_1D-N000.txt of set A' ...        
    'Log-Normal.mean0.50-sigma1.00_1D-N000.txt of set A' ...
    'Log-Poisson.h-0.50-coD1.00_1D-N000.txt of set B' ...
    'Log-Poisson.h-0.25-coD0.75_1D-N000.txt of set B' ...
    'Log-Normal.mean0.50-sigma1.00_1D-N000.txt of set B' ...
    'Log-Normal.mean0.25-sigma0.33_1D-N000.txt of set B' ...
    'Log-Poisson.h-0.50-coD1.00_1D-N000.txt of set C'...
    'Log-Normal.mean0.50-sigma1.00_1D-N000.txt of set C'...
    };

if ~exist('choice','var'),
    choice=1;
end;

if choice==1,        name='LogPoisson.A1.txt';
elseif choice==2,    name='LogNormal.A3.txt';
elseif choice==3,    name='LogPoisson.B1.txt';
elseif choice==4,    name='LogPoisson.B2.txt';
elseif choice==5,    name='LogNormal.B3.txt';
elseif choice==6,    name='LogNormal.B4.txt';
elseif choice==7,    name='LogPoisson.C1.txt';
elseif choice==8,    name='LogNormal.C3.txt';
elseif choice==0,    name='devilStaircase.txt';
end

fprintf('\nProcessing the file %s...', Nequiv{choice+1});

fid=fopen(name);
X=fscanf(fid,'%f');
fclose (fid);
N=length(X);

figure; plot(X);  
title('Original Signal');

Q = linspace(-5,5,51);
nvoice=10;
scale0=2; cut0=scale0;
noct=1;
noctave = floor(log2(N))-noct;
wavelet='Sombrero';

X(N+1:2*N) = X(N:-1:1);
noctave=noctave+1;
nscale  = nvoice * noctave;

% Computes the Continuous Wavelet Transform of the signal X
wt = CWT( X, nvoice, wavelet, noct, scale0 );

% Resize the matrix of WT coefficients
[n,nscale] = size(wt);
wt = wt(:, cut0*nvoice + 1 : nscale - nvoice);
[n,nscale] = size(wt);
scale = 2.^(log2(scale0) + (noctave-cut0-1)*(1:nscale)/nscale);
    
% Build Maxima Map relating the positions of the WTMM
par=100000000;
maxmap = MM_RWT( wt, par );
   
% Build Thermodynamic Partition Function
% Method 1:
Z = CalcThermoPartition( wt, maxmap, Q );
% Method 2: much longer !!!
if choice>=3, % used for the short series only
    Z2 = myCalcThermoPartition( wt, maxmap, Q );
else
    Z2=Z;
end

loscale=scale(1);
hiscale=scale(nscale);
% Method 1:
tau = FracScalExp( Z, scale, loscale, hiscale );
tau = -tau;  
% Method 2:
tau2 = myCalcMomentGenFun( Z2, scale, loscale, hiscale );
    
% Method 1:
HMIN=-2.; HMAX=1.5;
h = linspace( HMIN, HMAX, 101 );
dh = FracSingSpect( tau, Q', h );
i = find(dh >= 0);
dh = dh(i);    h = h(i);
% Method 2:
[ h2, dh2 ] = myCalcSingSpectrum ( tau2, Q );

figure; plot( Q, tau,'b+' ); hold on, plot( Q, tau2, 'r-' ); 
title('Multifractal Exponents');
xlabel('q');	ylabel('Tau(q)');
 
figure; plot( h, dh,'b+' ); hold on, plot( h2, dh2, 'r-' );
title('Multifractal Spectrum');
xlabel('h');	ylabel('D(h)');

fprintf('\n i) blue ''+'': results obtained processing the WTMM analysis scheme with the');
fprintf('\n    standard functions of the WaveLab library')
fprintf('\n ii) red ''-'': results obtained processing the WTMM scheme with functions:')
if choice>=3,
    fprintf('\n     - myCalcThermoPartition for the computation of the partition function');
end
fprintf('\n     - myCalcMomentGenFun for the computation of the exponents tau(q)');
fprintf('\n     - myCalcSingSpectrum for the computation of the spectrum (h,d(h))');