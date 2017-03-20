function [phi, y, s] = GenSparseChirps4(t, M, N, S)
% GenSparseChirps4.m
% [phi, y, s] = GenSparseVectors4(M, N, S)
%
% For testing of sparse algorithms.
% This function produces S non-zero coefficients in an N dimensional
% coefficient space from M measurements in 4 dimensions.
%
% Function output is an N X 4 matrix s, with S randomly assigned
% non-zero entries. It also produces a random sensing array, phi,
% which is created from M measurements of phi*s.
%
% Last modified: Ra Inta, Nov 15, 2012
%

% Initialize random number generator
    %randn('state',0); rand('state',2); % NOTE: This is the old syntax!!!
    rng('default')
    % Create sparse vector s and observation vector y
    
    s  = zeros(N, 4);
    
    for i = 1:4;
      sIdx = randperm(N);
      injections(i, :) = sIdx(1:S);       % Randomise S positions within s
    end
    
    
    f0 = 100 + injections(1, :);  % Do some scaling to get realistic injection frequency etc.
    tau0 = min(t) + (( max(t) - min(t) )/N)*injections(2, :);
    Q0 = (1000/N)*injections(3, :);
    d0 = (-20000/N)*injections(4, :);
    
    y = zeros(size(t));
    for i = 1:S;
        y = y + chirpxform(t, f0(i), tau0(i), Q0(i), d0(i));  % Exploit superposition property
    end
    
    phi = matrix_normalizer_array(randn(M, N, 4)); % Make phi completely random
    
    for i=1:4;
        Y(:,i) = phi(:,:,i)*s(:,i);
    end