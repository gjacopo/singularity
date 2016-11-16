% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% MYCONTWTSPEC - Estimates the multifractal Legendre spectrum of a 1-D 
% signal from the wavelet coefficients of a L2 continuous decomposition
% This is inpired by the function CONTWTSPEC of the library FracLab
% See function MYCONTWTPART
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [ h, dh, tau ] = mycontwtspec(part,scale,Q,ChooseReg,varargin)
% Usage
%     [ h, dh, tau ] = mycontwtspec(part,scale,Q,ChooseReg)
% Input parameters
%   o  part : real matrix [N_scale,N_Q], partition function
%   o  scale : real vector  [1,N_scale]. Analyzed scale vector
%   o  Q :  real vector [1,N_Q]. Exponents of the partition function
%   o  ChooseReg : 0/1 flag or integer vector [1,N_reg], (N_reg <=
%      N_scale) 
%         ChooseReg = 0 : full scale range regression
%         ChooseReg = [n1 ... nN_reg] : scale indices for the linear
%         regression of the partition function.
% Output parameters
%   o  h : Real vector [1,N_alpha], N_alpha <= N_Q
%      Singularity support of the multifractal Legendre spectrum
%   o  dh : real vector [1,N_alpha]. Multifractal Legendre spectrum
%   o  tau : real vector [1,N_Q]. Multifractal Exponents

nscale = length(scale); 
logscale = log2(scale(:)) ;

logpart = log2(part);

if isempty(varargin)
   varargin{1} = 'ls' ;
end
 
if ChooseReg == 1 % interactive mode
  
  figure('Tag','graph_reg')
  ireg = 1:nscale ;
  reg_log = ireg ;
  sth=findobj('Tag','StaticText_error');
  set(sth,'String','Press return to end !');
  
  while ~isempty(reg_log)   
    
      subplot(121); plot(logscale,logpart,logscale,logpart,'+'), axis tight
      ylabel('Partition function log_2(S_n(q))'), xlabel('log_{2}(scale)')

      reg_log = fracginput(2) ;
    
      if ~isempty(reg_log)
          ireg = find(min(reg_log(1,1),reg_log(2,1)) <= logscale & ...
          logscale <= max(reg_log(1,1),reg_log(2,1))) ;
      end
        
   
      for nq = 1:length(Q) 
          switch varargin{1}
              case 'ls'
                  slope = polyfit(logscale(ireg),logpart(ireg,nq),1) ;
              case 'wls'  
                  RegWeight = varargin{2}(ireg)./sum(varargin{2}(ireg)) ;
                  slope = monolr(logscale(ireg),logpart(ireg,nq),varargin{1},RegWeight) ;
              case {'pls'}
                  slope = monolr(logscale(ireg),logpart(ireg,nq),varargin{:}) ;
              case {'ml','lapls'}
                  slope = monolr(logscale(ireg),logpart(ireg,nq),varargin{:}) ;
              case {'linf','lsup'}
                  slope = newnewreg(logscale(ireg),logpart(ireg,nq),varargin{:},0) ;   
          end
  
          tau(nq) = slope(1) - 1 - Q(nq)/2 ;
      end
 
      % computation of the Legendre spectrum
      [dh,h] = flt(Q,tau) ;
      subplot(222), plot(h,dh) ; title('spectrum') ;
      xlabel('h'),ylabel('D(h)')
      subplot(224), plot(Q,tau) : title('\tau(q)') ; xlabel('q') ;
      
  end
       
  set(sth,'String','');

elseif ChooseReg == 0
  ireg = 1:nscale ;
  
elseif length(ChooseReg) > 1
  ireg = ChooseReg;
end

for nq = 1:length(Q) 
  
  switch varargin{1}
    case 'ls'
      slope = polyfit(logscale(ireg),logpart(ireg,nq),1) ;
    case 'wls'  
      RegWeight = varargin{2}(ireg)./sum(varargin{2}(ireg)) ;
      slope = monolr(logscale(ireg),logpart(ireg,nq),varargin{1},RegWeight) ;
    case {'pls'}
      slope = monolr(logscale(ireg),logpart(ireg,nq),varargin{:}) ;
    case {'ml','lapls'}
       slope = monolr(logscale(ireg),logpart(ireg,nq),varargin{:}) ;
    case {'linf','lsup'}
	    slope = newnewreg(logscale(ireg),logpart(ireg,nq),varargin{:},0) ;   
  end
  
  tau(nq) = slope(1) - 1 - Q(nq)/2 ;
  
end

% computation of the Legendre spectrum

[ dh, h ] = flt(Q,tau) ;



