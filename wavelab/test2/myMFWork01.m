% MFWork01: Multifractal Workout: CWT Analysis of Cantor Signal

	% modifiable parameters
	N = 1024;     % signal length; fairly large
	nvoice = 10;  % Following JS Bach, Well-Tempered Klavier

	% create a Brownian
	CantorMeasure = MakeFractal(N,3,'Deterministic',[.5 0 .5]);
	Devil  = cumsum(CantorMeasure); t = (.5:(N-.5))./N;
	figure; plot(t, Devil); title(sprintf(' Devil Staircase Signal'));
%fi=fopen('devil_staircase.txt','w');
% fprintf(fi,'%d %f \n',[1:N;Devil]);
%fprintf(fi,'%f\n',Devil);
%fclose(fi);

% make CWT
	% Devil = Devil - Devil(N) .*t;	
%	Devil_cwt = CWT(Devil,nvoice,'Sombrero');
	[Devil_cwt,lstscale] = myCWT(Devil,nvoice,'Sombrero');
	sz = size(Devil_cwt); nscale = sz(2);

	% display CWT
	figure; ImageCWT(Devil_cwt,'Individual','hot');
	title('CWT');
	
	% Build Maxima Map
	Devil_maxmap = myWTMM(Devil_cwt);
	% display maxmap
	figure; ImageWTMM(Devil_maxmap)
	
%	Devil_maxmap = WTMM_2(Devil_cwt);
%	% display maxmap
%	figure; ImageWTMM(Devil_maxmap);

    qlist = linspace(-5,5,51);
    
    z = myCalcPartitionFunc(Devil_cwt, Devil_maxmap, qlist);
    figure, imagesc(z), colormap jet
% z = CalcThermoPartition( Devil_cwt, Devil_maxmap, qlist );
% nscales=length(lstscale);
% z=z(:,10:nscales);
% lstscale=lstscale(10:nscales);
%fi=fopen('partfunc-wave.fi','w');
%fprintf(fi,'%g ',z);
%fclose(fi);

 figure, PlotThermoPartition( [], z, Devil_cwt, lstscale, qlist );
%save 'partfunc-wtmm.fi' z -ascii -double
%figure, imagesc(z), colormap gray

tau = myCalcMomentGenFun( z, lstscale);
 
%      z = CalcThermoPartition(Devil_cwt, Devil_maxmap, qlist);
  figure, PlotThermoPartition( [], z, Devil_cwt, lstscale, qlist );
  
%tau2 = CalcMomentGenFun( z, lstscale);
        figure, plot(qlist,tau,'b+');
        %hold on; plot(qlist,tau2,'y-');
       
       
 [ h, dh ] = myCalcSingSpectrum ( tau, qlist );
    figure; plot(h,dh,'b+');
 return 

 d = CalcGenFracDimen( z, qlist, lstscale );
 figure; PlotGenFracDimen( d, qlist,lstscale  );


f = CalcFracSpectrum( tau, qlist' );
figure; PlotMultiSpectrum( f );

