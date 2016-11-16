% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% MYCALCMOMENTGENFUN -- Calculate Moment Generating Function
%  Usage
%    tau = myCalcMomentGenFun(z,scale,loscale,hiscale)
%  Inputs
%    z         matrix nexp by nscale of z(q,a) ``Thermo Partition Func''
%    scale         list of scales 
%    loscale   optional min scale to fit line by
%    hiscale   optional max scale to fit line by 
%  Outputs
%    tau       vector 1 by nscale of moments
%
%  Description
%    tau(q) = Slope [ log(z(q,a))  versus log(a) ]
%
%  See Also
%    CalcThermoPartition
%
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function tau = myCalcMomentGenFun(z,scale,loscale,hiscale)

	if nargin < 3,
		loscale = 0.;
		hiscale = 10^32;
	end
	
	[nq,nscale]     = size(z);
	tau = zeros(nq,1);
	if nscale ~= size(scale),
		disp('CalcMomentGenFun: no match between z and a')
	end
	
    % to be sure that scale is a line vector 
   if (size(scale,1)>1) scale = scale'; end
   
	window = (loscale <= scale) & (scale <= hiscale);
	ix     = find(window);
    % number of scales used
    s0 = sum(window);
    
 	if length(ix)>1,
   
       lz =  log2(z(:,ix))';
       lsc = log2(scale(ix));
      
        s1 = sum( lsc );
        s2 = sum( lsc .* lsc );
  
        t0 = sum( lz );
        t1 = lsc * lz;
  
        tau = (t1*s0-s1*t0) / (s0*s2-s1*s1);
    
	else
		disp('CalcMomentGenFun: not enough data to determine slope')
	end 	 
    
    % this is one  trick of the wavelab processing:
    tau = -tau;
