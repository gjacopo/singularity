% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% myCalcThermoPartition -- Build Thermodynamic Partition Function
%  Usage
%    z = myCalcThermoPartition(wt,mm,Q)
%  Inputs
%    wt        continuous wavelet transform output by CWT
%    mm        maxmap output by MMWT
%    Q         list of moments
%  Outputs
%    z         matrix nexp by nscale of z(q,a)
%  Description
%    z(q,a) = sum_i( |CWT(a,b(i))|^q ),  where b = (b(i)) is a list
%    of wavelet transform maxima  
%
% Uses the function BuildSkelMapFast and PruneSkelMap to refine the skeleton
% of the WTMM representation.
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function  z = myCalcThermoPartition(wt,mm,q)
  if nargin < 6,
	oct = 2;
	sc = 4;
  end
  
  nq = length(q);
  if size(q,2)>1 q=q'; end;
  
    [n,nscale] = size(wt);
	noctave = log2(n)-oct;
	nvoice  = nscale/noctave;

    % Identify ridges
	[skm,skp,skl] = BuildSkelMapFast( mm );

	% Eventually prune ridges
	[skm,skp,skl] = PruneSkelMap(wt,.001,1,skm,skp,skl);
    % Note: in fact, we only need the list of (scale,position) of the extrema
    
    nmax = length(skm);
    nmax = nmax / 2.;
    
 	vec = zeros(2,nmax);
	% Transform the vector skellist in matrix (npunt,2) 
    vec(:) = skm;
    vec = vec';
    % Now:  
    %    - vec(:,1) is the list of the scales where the maxima can be 
    %      found, disregarding of the chain they belong to
    %    - vec(:,2) is the list of the position of the maxima values at
    %      the corresponding scale, also disregarding of the chain they 
    %      belong to
    % So, note that we don't make any distinction between the different 
    % chain of appartenance of the maxima (we dont care anymore about the 
    % chain the maxima value is computed on, we only care about the scale
    % it was computed at)

    % Compute the absolute value of the CWT over the different extracted
    % maxima 
    sc = vec(:,1);
    ipos = vec(:,2);
    wt = ShapeAsRow(wt)';
    ind = (sc - 1) * n  + ipos;
    W = [ sc, wt(ind) ];
    
    z=zeros(nq,nscale);
    for k=1:nscale,
	    j = find(W(:,1)==k);
		if ~isempty(j),
			for i=1:nq,
				z(i,k) = sum(abs(W(j,2)) .^ q(i));
            end     
        else
            z(:,k) = eps.^q;
		end

    end
                
  


