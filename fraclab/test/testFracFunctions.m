% TESTFRACFUNCTIONS - Performs the multifractal analysis of some determined
% test signals with the modified FracLab codes.
% Usage
%       testFracFunctions(choice)
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
function testFracFunctions(choice)

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
scale0=2;
noct=1;
noctave = floor(log2(N))-noct;

X(N+1:2*N) = X(N:-1:1);
noctave=noctave+1;
nscale  = nvoice * noctave;

nscale  = nvoice .* noctave;

fmin = 2^(-noctave);
fmax = min( 0.5, 2^(-log2(scale0)));

% Mode='L2'
[wt,scale] = contwt(X,fmin,fmax,nscale,0) ; 
% Mode='L1'
% [wt,scale] = cwt1d(X,fmin,fmax,nscale,0) ; 
% !!! Note : other bug FracLab - the function cwt1d crashes with
% large time-series. !!!
         
cut0=scale0;
wt = wt(cut0*nvoice+1 : nscale - nvoice,:);       
scale = scale(cut0*nvoice+1 : nscale - nvoice);        
[n,nscale] = size(wt);

% Use the full range of scales for the estimation
% Mode='L2'
part = mycontwtpart( wt, scale, Q, 1);
[ h, d, tau ] = mycontwtspec( part, scale, Q, 0 );
% Mode='L1'
% part = mycwtpart( wt, scale, Q, 1);
% [ h, d, tau ] = mycwtspec( part, scale, Q, 0 );

% display
figure; subplot 121, plot( Q, tau );
title('Exponents - Full range');
xlabel('q');	ylabel('Tau(q)');
subplot 122, plot( h, d );
title('Spectrum - Full range');
xlabel('h');	ylabel('D(h)');
            
% Make the estimation for different parameters
fprintf('\nPress return to finish');
% Mode='L2'
[ h, d, tau ] = mycontwtspec( part, scale, Q, 1 ) ;
% Mode='L1'
% [ h, d, tau ] = mycwtspec(wt,scale,Q,1,'ls') ;
