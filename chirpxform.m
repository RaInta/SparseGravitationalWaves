function psi = chirpxform(t, f, tau, Q, d);


% ChirpXform.m
%
% psi = chirpxform(t, f, tau, Q, d);
% This performs the chirp xform, as described in Chassande Mottin et al
% CQG 27 194017 (2010)
% I would like to calculate the mutual coherence between the chirplet and the
% time basis.
% 
% Created: Nov 4, 2011 Ra Inta
% Last modified: Jan 18, 2012 R.I.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to replicate the example in the paper, take these values
%
% f = 350;   % Centre frequency in Hz
% tau = 0.05; % Centre time
% Q = 50;    % Quality factor
% d = -5000; % Modulation rate, Hz/s
% t = 0:0.001:0.1;  % Time in seconds
%
% If you give this function only one argument, these values are adopted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 1
    f = 350;   % Centre frequency in Hz
    tau = 0.05; % Centre time
    Q = 50;    % Quality factor
    d = -5000; % Modulation rate, Hz/s
    t = 0:0.001:0.1;  % Time in seconds
end

A = (8*pi*f^2/Q^2)^0.25;  % Normalisation factor, so that \int(norm(psi)^2) = 1

psi = A*exp(-2*pi*(f/Q)^2 * (tau - t).^2).*exp(2*pi*1i*(f*(tau - t) + 0.5*d*(tau - t).^2 ));

return


