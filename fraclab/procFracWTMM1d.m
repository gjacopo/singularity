% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% PROCFRAWTMM1d -- Script that launches the analysis of multifractal signals
% Usage:
%                 procFracWTMM1d( set, process, Q, R ) 
% Input parameters:
%     - set : set of data to be processed: 'A', 'B', 'C' or 'all' 
%     - process :  integer flag in {0,...,4}:
%           process=0 or 'all' : all processed
%           process=1 : Log-Poisson.h-0.50-coD1.00_1D only
%           process=2 : Log-Poisson.h-0.25-coD0.75_1D ____
%           process=3 : Log-Normal.mean0.50-sigma1.00_1D only
%           process=4 : Log-Normal.mean0.25-sigma0.33_1D ____
%     - Q : list of moments (def.: the Q list of the article)
%     - R : ranges of scales to perform the tau(q) estimation (def.: the
%       range [4,100] of the article)
% See also MYANALYZEFRACSERIESWTMM
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Note regarding the installation:
% Verify first that the FracLab library has benn installed
% Two possibilities:
% 1) Use the functions of the subdirectory fraclib/ of the current directory.
% 2) Add the path of the functions of the Fraclab libraries to the Matlab PATH 
% variable. This can be done thanks to the command ADDPATH, e.g. typing 
% something like: 
%   > addpath(genpath('C:\MyPathToFraclab\Fraclab'));
% where in my case C:\MyPathToFraclab\Fraclab is the path of the directory 
% containing the Fraclab functions.
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function procFracWTMM1d( set, process, Q, R ) 

% check parameters
if ~exist('set','var'), 
    set='all'; 
elseif isempty(strmatch(set,strvcat('A','B','C','all'),'exact')),
    error('Parameter ''set'' must be ''A'', ''B'', ''C'' or ''all'' (see help)');
end;
if ~exist('process','var'), 
    process=0; 
elseif strcmp(num2str(process),'all'), 
       process=0; 
elseif process<0 | process>4 
       error('Parameter ''process'' must be 0, 1, 2, 3 or 4 (see help)'); 
end;

%% Parameters from the article
% List of moments q for computing the partition function Z(q,a)
if ~exist('Q','var'), 
    Q = [ -4. -3.6 -3.2 -3. -2.8 -2.6 -2.4 -2.2 -2. ...
        -1.8 -1.6 -1.4 -1.2 -1.1 -1. ...
        -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0. ...
        0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 ...
        0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1. 1.05 ...
        1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2. ...
        2.3 2.6 3. 3.5 4. 5. 6. 7. 8. ];
end;
% Range of regression on scales a for etimating the exponents tau(q)
if ~exist('R','var'),
    R = [4,100];
end;
% Make tests for R=0: the results are probably better.

%% Some parameters for performing the WTMM
% Number of octaves discarded by the analysis
noct=1;
% The number of octaves for computing the CWT will be:
%    noctave = floor(log2(N))-noct;
% with N the lenght of the serie. Then I assign the value:
%    fmin = 2^(-noctaves)
% for the minimum frequency of analysis. If you want, you can eventually 
% change this parameter in the function myAnalyzeSeriesFracWTMM.

% Minimal scale of analysis
scale0=1.5;

% Number of voices per octave
nvoice=10;
% The total number of scales of analysis will be
%    nscales  = nvoice .* noctave;

% Mode: important parameter which defines which CWT to be used
mode='L2';
% !!! It seems that there is a problem with the function cwt1d
% for large series (only the library is provided) !!!
% Thus keep the mode='L2' to use functions: CONTWT, MYCONTWTPART
% and MYCONTWTSPEC.

%% Some dummy but useful variables
benchmark='benchmark';

LDIR = {'A' 'B' 'C'};
DIM = [16384 4096 1024];
NSERIES = [1 10 100];
if ~strcmp(set,'all'),  % process one set only
    ianal = strmatch(set,strvcat('A','B','C'),'exact');
    nanal=ianal;
else
    ianal=1;
    nanal=3;
end

LF = {'Log-Poisson.h-0.50-coD1.00_1D' 'Log-Poisson.h-0.25-coD0.75_1D' ...
    'Log-Normal.mean0.50-sigma1.00_1D' 'Log-Normal.mean0.25-sigma0.33_1D'};
if process~=0,   % process one type of multifractal process only
    iproc=process;
    nproc=iproc;
else
    iproc=1;
    nproc=4;
end;

% Some parameters for function myAnalyzeFracSeriesWTMM
binary = 0; % if text files, otherwise put 1

for i=ianal:nanal,
    % for each set: 1x16384, 10x4096 or 100x1024
    
    % enter the name of the input/output directories
    d = LDIR{i};
    storedir = sprintf('../%s/%s', benchmark, d);
    savedir = sprintf('%s/', d);
    
    % properties of the analyzed dataset
    N = DIM(i);
    nseries = NSERIES(i);
  
    for j=iproc:nproc,
        % for each multifractal process: Log-Poisson or Log-Normal

        % Generic name of the files under study    
        base = LF{j};

        % Main computation
        [ h, d, tau ] = myAnalyzeSeriesFracWTMM( storedir, savedir, base, binary, ...
                nseries, N, Q, R, scale0, nvoice, noct, mode);    
            
    end 

end
