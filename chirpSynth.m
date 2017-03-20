% chirpSynth.m
%
% This script creates a time series to simulate a gravitational wave chirp
% for the purposes of compressive sampling analysis, in an attempt to obtain 
% superresolution in position reconstruction
%
% Created: October 27, 2011 Ra Inta
% Last modified: October 27, 2011 R.I.

fs = 2^14; % 16 kHz sampling rate 
t = 0:1/fs:10; % 10 second data

N = length(t);

% Construct chirp time series
y = 1e-5*ones(1,N);
t0 = 3*fs; 
t1 = 8*fs;
y(t0:t1) = chirp(t(t0:t1), 100, 8, 3000);
n = randn(1,N);
n = abs(n./max(abs(n)));  % Normalise noise and make positive
y = y + 1e-5*n; % Additive Gaussian noise

specgram(y,256,fs,256,250);  % Plot the spectrogram
colorbar

% Find 3rd order Renyi entropy (requires S_renyi Matlab function)

disp(sprintf('\n>> Renyi entropy of time series: %0.3f', S_renyi(y, 3)))
disp(sprintf('>> Renyi entropy of background noise: %0.3f \n', S_renyi(n, 3)))