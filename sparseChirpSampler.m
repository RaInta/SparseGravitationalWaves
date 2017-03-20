function sparseChirp = sparseChirpSampler(t, f0, tau0, Q0, d0);

% sparseChirpSampler.m
%
% sparseChirp = sparseChirpSampler(t, f0, tau0, Q0, d0);
% 
% This is meant to sample randomly from central parameters f0, tau0, Q0 
% and d0 from the chirpxform.m function.
%
%
% Created: Aug 28, 2012 Ra Inta
% Last modified: Aug 28, 2012 R.I.
%
%

% Old version:
%function sparseChirp = sparseChirpSampler(t, f0, tau0, Q0, d0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These should eventually be input parameters in the function
fRange = 700;
tauRange = 0.1;
QRange = 100;
dRange = 10000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test values:
% f0 = 350;   % Centre frequency in Hz
% tau0 = 0.05; % Centre time
% Q0 = 50;    % Quality factor
% d0 = -5000; % Modulation rate, Hz/s
% t = 0:0.001:0.1;  % Time in seconds
% Pastable version:
% f0 = 350; tau0 = 0.05; Q0 = 50; d0 = -5000; t = 0:0.001:0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In case you just want a quick example
if nargin == 0
    disp('Warning - using default values for demonstration purposes')
    f0 = 350;   % Centre frequency in Hz
    tau0 = 0.05; % Centre time
    Q0 = 50;    % Quality factor
    d0 = -5000; % Modulation rate, Hz/s
    t = 0:0.001:0.1;  % Time in seconds
end

fMin    = max(0, f0 - fRange/2);
tauMin  = max(0, tau0 - tauRange/2);
QMin    = max(0, Q0 - QRange/2);
dMin    = max(0, d0 - dRange/2); 

f   = fMin + fRange*rand(1,1);
tau = tauMin + tauRange*rand(1,1);
Q   = QMin + QRange*rand(1,1);
d   = -( dMin + dRange*rand(1,1) ); % Note that d < 0


sparseChirp = chirpxform(t, f, tau, Q, d);

return
