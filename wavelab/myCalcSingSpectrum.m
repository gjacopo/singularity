% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% MYCALCSINGSPECTRUM -- Compute the spectrum (h,D(h)) from the multifractal
% exponents tau(q)
% Usage
%            [ h, dh ] = CalcSingSpectrum ( tau, Q )
% Input
%     tau    multifractal exponents output by CalcMomentGenFun
%     Q      list of moments used for estimation  
% Output
%     h      singularity exponents' values
%     dh     singularity spectrum
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [ h, dh ] = myCalcSingSpectrum ( tau, Q )
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if ~exist('dec','var'),
    dec=0;
end;

nq = length(Q);
   
tau = ShapeAsRow( tau );
Q = ShapeAsRow( Q );
% 
tau = tau - Q/2;

% First compute the singularity exponent h:
%        h = \frac{d\tau}{dq}

ltau=lshift(tau);  rtau=rshift(tau);
qlist = ShapeAsRow( Q );
lq=lshift(Q);  rq=rshift(Q);

h = (ltau - rtau) ./ (lq - rq);
h = h(2:nq-1);
% note: this is equivalent to :
% for iq=2:nq-1
%    h[iq] = (tau[iq+1]-tau[iq-1])/(qlist[iq+1]-qlist[iq-1]);
% end


% Then compute the density spectrum:    
%        D(h) = q h -tau(q) 
dh = qlist(2:nq-1) .* h - tau(2:nq-1);  
 
