% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% MYCWTPART - Estimates the partition function of a 1-D signal from the 
% wavelet coefficients of a L2 continuous decomposition.
% This is inpired by the function CWTSPEC of the library FracLab.
% The modification regards the parameter part storing the partition function:
%    i) the function itself is returned, not his logarithm
%    ii) it is returned as the first output
% See function MYCWTSPEC.
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [part, maxmap]  = mycwtpart( wt, scale, Q, FindMax )
% Usage
%        [part, maxmap] = mycwtpart( wt, scale, Q, FindMax )
% Input parameters
%   o  wt : Real or complex matrix [N_scale,N]
%      Wavelet coefficients of a continuous wavelet transform (output of
%      contwt or contwtmir))
%   o  scale : real vector  [1,N_scale]. Analyzed scale vector
%   o  Q :  real vector [1,N_Q]. Exponents of the partition function
%   o  FindMax : 0/1 flag.
%         FindMax = 0 : estimates the Legendre spectrum from all coefficients
%         FindMax = 1 : estimates the Legendre spectrum from the local Maxima
%         coefficients of the wavelet transform
%         Default value is FindMax = 1
% Output parameters
%   o  part : real matrix [N_scale,N_Q], partition function
%   o  maxmap : binary matrix of position of WTMM

nscale = length(scale) ;

if FindMax == 1 
  maxmap = findWTLM(wt,scale) ;
elseif FindMax == 0
  maxmap = ones(size(wt)) ;
end

% matrix reshape of the wavelet coefficients

detail = abs(wt.') ; 
maxmap = maxmap.' ;

% computation of the partition function and the mass function

for nq = 1:length(Q)
  for j=1:nscale
    max_idx = find(maxmap(:,j) == 1) ;
    DetPowQ = detail(max_idx,j).^Q(nq) ;
    part(j,nq) = mean(DetPowQ);
  end
end
