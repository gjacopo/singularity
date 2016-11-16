
[n,b] = size(sig); % n = no. of datapoints 
J = log2(n); % J = total no of dyadic levels. 2^J=n 
alpha = FWT_ATrou(sig,0); % perform standard Algorithm A Trous Wavelet Transform 
u = MM_DWT(alpha,0); % u = modulus maxima identified with 1 
a = u.*alpha; % a = modulus maxima identified with maxima value 

