% test
function scale = makeyourtest()
choice=1;
if choice==1
    % signal 1
  %  fid=fopen('Log-Poisson.h-0.50-coD1.00_1D-N000.txt','r');
  fid=fopen('Log-Poisson.h-0.50-coD1.00_1D-N000.txt','r');
    X=fscanf(fid,'%f');
    fclose (fid);
    N=length(X);
    
elseif choice==2
    % signal 2
    N=4096;
    CantorMeasure = MakeFractal(N,3,'Deterministic',[.5 0 .5]); 
    X = cumsum(CantorMeasure); t = (.5:(N-.5))./N;
    
elseif choice==3
    %signal 3
    fid=fopen('devil_staircase.txt','r');
    X=fscanf(fid,'%f');
    fclose (fid);
    N=length(X);
end


Q = linspace(-10,10,61);
loscale=4; hiscale=100;
nvoice=12;
scale0=4;
cut0=4;
noct=2;
%wavelet='Sombrero';
wavelet='DerGauss';

  %  fmin = 2^(-noctave)
    % fmin = 2^(-9);
   % fmax = min( 0.5, 2^(-log2(scale0)))
    % [ wt, scale ] = contwtmir( X, fmin, fmax, nscale, 0 );
       
%%   

% Miror the signal
X(N+1:2*N) = reverse(X);
	
figure; plot(X);  
title('Original Signal');

wt = RWT(X,nvoice,wavelet);	
wt2 = CWT( X, nvoice, wavelet );
  noctave = floor(log2(N))-noct +1;
nscale  = nvoice .* noctave;

[n,nscale] = size(wt);
 wt = wt(:, cut0*nvoice+1 : nscale - nvoice);
 [n,nscale] = size(wt);
 scales = 2.^(log2(scale0)+(noctave-cut0+1)*(1:nscale)/nscale);
scales2=scales;

  par=100000000;
maxmap = MM_RWT(wt,par);
z = FracPartition(wt,maxmap,Q);   

maxmap2 = MM_RWT(wt2,par);

z2 = z;
%myCalcThermoPartition( wt2, maxmap2, Q );

PlotThermoPartition(linspace(min(Q),max(Q),7),z,wt,scales,Q);
hold on, PlotThermoPartition(linspace(min(Q),max(Q),7),z2,wt2,scales2,Q);

% comparer
tau = FracScalExp(z,scales);
tau = -tau;
tau2 = myCalcMomentGenFun( z2, scales2 );
figure, PlotMomentGenFun(tau,Q)
hold on, plot( Q, tau2,'r-' );

alpha = linspace(-1.,.9,101);
f = FracSingSpect(tau,Q',alpha);
i=find(f>=0);
f=f(i);
alpha=alpha(i);
figure,	
 PlotMultiSpectrum(f,alpha);
 	hold on
 % axis([.4 .9 0 0.7]);

  [ h, dh ] = myCalcSingSpectrum ( tau, Q );
  plot( h, dh,'r-' );	 
  
hold off   
   
return;

 %% METHODE WAVELAB
  [wt,scale] = myCWT( X, nvoice, wavelet );
  sz = size( wt ); nscale = sz( 2 );
 
  scale2 = 2.^(log2(scale0)+noctave*(1:nscale)/nscale);
 
 figure, plot(scale); hold on, plot(scale2);
 hold off
 
 maxmap = WTMM( wt );
 par=100000000;
 maxmap2 = MM_RWT(wt,par);        
 figure; ImageWTMM( fliplr(maxmap) );
 figure; ImageWTMM( fliplr(maxmap2) );

 z = myCalcThermoPartition( wt, maxmap, Q );
 tau = myCalcMomentGenFun( z, scale, loscale, hiscale );
[ h, dh ] = myCalcSingSpectrum ( tau, Q );

 figure; plot( Q, tau );
 xlabel('q');	ylabel('Tau(q)');	
 title('Multifractal Exponents - 1'); 
 figure; plot( h, dh );	 
 title('Multifractal Spectrum - 1'); 
 xlabel('h');	ylabel('D(h)');

  z = myCalcThermoPartition( wt, maxmap2, Q );
 tau = myCalcMomentGenFun( z, scale2, loscale, hiscale );
[ h, dh ] = myCalcSingSpectrum ( tau, Q );
 figure; plot( Q, tau );
 xlabel('q');	ylabel('Tau(q)');	
 title('Multifractal Exponents - 2'); 
 figure; plot( h, dh );	 
 title('Multifractal Spectrum - 2'); 
 xlabel('h');	ylabel('D(h)');


