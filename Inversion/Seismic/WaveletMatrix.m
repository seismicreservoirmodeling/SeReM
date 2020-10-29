function W = WaveletMatrix(wavelet, nsamples, ntheta) 

% WAVELET MATRIX computes the wavelet matrix for discrete convolution
% INPUT w = wavelet
%       nsamples = numbr of samples
%       ntheta = number of angles
% OUTUPT W = wavelet matrix

% Written by Dario Grana (August 2020)

W = zeros(ntheta*(nsamples-1));
[~, indmaxwav] = max(wavelet);  
for i=1:ntheta
    wsub = convmtx(wavelet', (nsamples-1))';
    indsub = (i-1)*(nsamples-1)+1:i*(nsamples-1);
    W(indsub, indsub) = wsub(indmaxwav:indmaxwav+(nsamples-1)-1,:);
end