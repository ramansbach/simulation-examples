function a = akaike(k,n,r)
%function that computes the akaike information criterion
%k is the number of parameters, n is the number of data points, r is the
%sum squared error
a = 2*k + n*log(r/n);