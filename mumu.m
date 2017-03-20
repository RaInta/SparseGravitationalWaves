unction [muCoherence] = mumu(phi, psi)

% [muCoherence] = mumu(phi, psi)
%
% Mutual coherence, used in compressive sampling methods to establish
% the number of measurements required to take in order to guarantee
% accurate reconstruction of sparse signals:
%
% M > C mumu.^2 S log(N)
% where C is a constant, S is the underlying sparsity of the model (signal)
% vector and N is the dimensionality of the original measurement space.
%
% Created: October 29, 2011 Ra Inta
% Last modified: Nov 18, 2011 R.I.


[M N] = size(phi);
muCount = zeros(1, M);

U = phi*psi; 
for iIdx = 1:N
    for jIdx = 1:N;
        if iIdx ~= jIdx;
            muCount(iIdx, jIdx) = max( max( abs( U(:, iIdx) .* U(:, jIdx) )./( norm( U(:, iIdx) ) * norm( U(:, jIdx) ) )) );
        else
            muCount(jIdx, jIdx) = 0;
        end
    end
end
muCoherence = max(max(muCount));

