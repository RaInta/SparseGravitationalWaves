% chirpxform_variables.m
%
%
% This loads variables so that the Chirp transform function, chirpxform.m
% , makes sense, from Chassande Mottin's paper on the chirplet transform
% in GW data analysis.
% 
% Created: Nov 4, 2011, Ra Inta
% Last modified: Nov 4, 2011, R.I.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to replicate the example in the paper, take these values

% f = 350;   % Centre frequency in Hz
% tau = 0.05; % Centre time
% Q = 50;    % Quality factor
% d = -5000; % Modulation rate, Hz/s
% t = 0:0.001:0.1;  % Time in seconds


f = 100;   % Centre frequency in Hz
tau = 0.05; % Centre time
Q = 50;    % Quality factor
d = -10000; % Modulation rate, Hz/s
t = 0:0.001:0.1;  % Time in seconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


psi = chirpxform(t, f, tau, Q, d);