% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% MYAPPROXLEGENDRE -- Computes the multifractal exponents
% Direct method to compute the spectrum thanks to the Legendre transform
% Used for comparison: this function is EQUIVALENT TO CALCMOMENTGENFUN
%     tau = ApproxLegendre (z,scale,loscale,hiscale)
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function tau = myApproxLegendre (z,scale,loscale,hiscale)
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% this is nothinig else that a simple 'translation' of the Legendre Transform
% algorithm implemented in the C program wtmm_processor
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if nargin < 3,
		loscale = 0.;
		hiscale = 10^32;
	end
	
	[nq,nscale]     = size(z);
	tau = zeros(nq,1);
	if nscale ~= size(scale),
		disp('CalcMomentGenFun: no match between z and a')
	end
	
    count = 0.;
    sumsc = 0.;
    sumsqsc = 0.;
    sumpf = zeros(nq,1);
    sumscpf = zeros(nq,1);
    
    for iws=2:nscale    
        % log of the current scale 

        if(scale(iws)>=loscale && scale(iws)<=hiscale) 
            Lws = log( scale(iws) );
            % number of scales considered till now 
            count = count + 1;
            % sum of scales 
            sumsc = sumsc + Lws;
            % sum of squared scales 
            sumsqsc = sumsqsc + Lws*Lws;
      
            for iq=1:nq 
	            % log of the partition function 
            	LZ = log( z(iq,iws) );
	            % sum of the partition function 
	            sumpf(iq) = sumpf(iq) + LZ;
	            % sum of partition function weightened by the scale  
            	sumscpf(iq) = sumscpf(iq) + Lws * LZ;	
            end
        end
    end % end loop over the scales
             
    for iq=1:nq 
            tau(iq) = (count*sumscpf(iq) - sumsc*sumpf(iq)) / ...
                       (count*sumsqsc - sumsc*sumsc);
    end
    
