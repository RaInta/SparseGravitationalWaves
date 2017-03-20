% OMPerrAnalysis.m
%
% This script takes a range of residual noise values to set
% as a residual vector for Orthogonal Matching Pursuit 
% based compressive sampling. This level, epsilon, may be used 
% to define false alarm rates/upper bounds for Chi-squared distributions.
%
% Requires: GenSparseVectors.m (which requires matrix_normalizer.m ) and
% OMP.m
%
% Created: October 28, 2011 Ra Inta
% Last modified: October 28, 2011 R.I.

M = 5E2;        % Number of measurements
N = 1000;        % Number of original coefficients
S = 40;         % Expected sparsity of vector

[phi, y, s] = GenSparseVectors(M, N, S);

for loopIdx = 1:7
    epsilon = 10^(loopIdx - 6);
    [sOMP, R] = OMP(phi, y, epsilon);
    Rmax(loopIdx) = max(abs(R));  % Record maximum error
    Rmin(loopIdx) = min(abs(R));  % Record minimum error
    sRecon{loopIdx} = sOMP; % Record reconstructed vector
    Rrecon{loopIdx} = R;
end
