function [phi, y, s] = GenSparseVectors(M, N, S)
% GenSparseVectors.m
% [phi, y, s] = GenSparseVectors(M, N, S)
%
% For testing of sparse algorithms.
% This function produces S non-zero coefficients in an N dimensional
% coefficient space from M measurements. 
% Function output is an N-length vector s, with S randomly assigned
% non-zero entries. It also produces a random sensing matrix phi,
% which is created from M measurements of phi*s.

% Initialize random number generator
    %randn('state',0); rand('state',2); % NOTE: This is the old syntax!!!
    rng('default')
    % Create sparse vector s and observation vector y
    sIdx = randperm(N);
    sIdx = sIdx(1:S);       % Randomise S positions within s
    s  = zeros(N, 1);
    s(sIdx) = randn(S, 1);  % Randomise amplitudes
    
    phi = matrix_normalizer(randn(M, N)); % Make phi completely random
    y = phi*s;