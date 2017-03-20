function [x, R] = OMP(phi, y, epsilon)

%  [x, R] = OMP(phi, y, epsilon)
%
% Orthogonal Matching Pursuit recursive reconstruction algorithm
% for compressive sampling.
%
% This solves the inverse problem y = phi*x, by finding x (N X 1 vector),
% where A = M X N measurement matrix, y = M X 1 measurement vector.
% epsilon is a scalar representing the stochastic/residual noise.
% This can be used to define a confidence limit to the Chi-squared distribution.
%
% Very similar to CLEAN routine used in radio interferometry; 
% the stopping criterion is when the residual vector R is reduced below the 
% threshold, epsilon
%
% Created: Ra Inta, October 13, 2011
% Last modified: Jan 18, 2012 R.I.



if nargin < 3
     epsilon = 1e-5;
     disp('>> Warning: setting residual error to default [1E-5]')
end

xTmp = zeros(size(phi, 2), 1);
r = y; L = []; Psi = [];  % Initialise vectors

while ( norm(r) > epsilon)
    l = phi'*r;   % Take pseudo-inverse of phi etc.
    [mostSig, mostSigIdx] = sort(abs(l));  % Select the most significant column
    sortBack = length(mostSig):-1:1;
    mostSig = mostSig(sortBack);   % Have to work around because this version doesn't have sort-descend
    mostSigIdx = mostSigIdx(sortBack);
    I = mostSigIdx(1);      % Take most significant column index
    L = [L' I']';           % Add to cumulative set of most significant indices 
    Psi = phi(:, L);        % Most significant column
    xTmp = Psi\y;           % Find reconstructed x from this column
    yApprox = Psi*xTmp;
    r = y - yApprox;
end

x(L) = xTmp;            % Final reconstructed vector
R = r;                  % Final residual vector