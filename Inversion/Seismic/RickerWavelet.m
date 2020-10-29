function [w,tw] = RickerWavelet(freq, dt, ntw)

% RICKER WAVELET computes the Ricker wavelet
% INPUT freq = dominant frequency
%       dt = time sampling rate
%       ntw = number of samples of the wavelet
% OUTUPT w = wavelet
%        tw = wavelet time

% Written by Dario Grana (August 2020)

tmin=-dt*round(ntw/2);
tw=tmin+dt*(0:ntw-1)';
w=(1-2.*(pi^2*freq^2)*tw.^2).*exp(-(pi^2*freq^2)*tw.^2);