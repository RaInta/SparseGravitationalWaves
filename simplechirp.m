function yChirp = simplechirp(t, w0, w1)

% yChirp = simplechirp(t, w0, w1)
%
% This function generates a simple chirp with base (angular) frequency
% w0 and chirp rate w1, applied to time vector t
%
%
% Created: Nov 28, 2011 Ra Inta
% Last modified: Nov 28, 2011 Ra Inta


yChirp = cos(2*pi*(w0*t + 0.5*w1*t.^2));