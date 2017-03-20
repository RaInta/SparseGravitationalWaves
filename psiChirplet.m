function yChirplet = psiChirplet(t, f0, tau0, Q, d) 

% yChirplet = psiChirplet(t, f0, tau0, Q, d)
%
% Make a normalised GW chirplet wave-form from time vector t
% with centre frequency f0, centre time tau0, Q-factor Q
% and chirp rate d
%
%%%%%%%%%%%%%%% Make chirplets %%%%%%%%%%%%%%%%%
% This is from Chassande-Mottin's Chirplet paper
%
% Example:  
%   dt=1e-3;
%   t=0:dt:0.3;
%   fs=1/dt
%   psi=psiChirplet(t, 100, 0.15, 500, 20);
%   specgram(psi,2^6,fs)
%   axis([0 0.3 0 200])
%
% Created: Nov 26, 2011 Ra Inta
% Last modified: Nov 26, 2011 R.I.

A = (8*pi*f0^2 / Q^2)^0.25; % Normalisation factor

Tmt = t - tau0;

yChirplet = A * exp(-(2*pi*f0)^2 * Tmt.^2 / Q^2) .* exp(2*pi*i*(f0*Tmt + 0.5*d*Tmt.^2));



