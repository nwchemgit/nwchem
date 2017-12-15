function [score] = eval_diff_ent_gau(tnsr)
    [x, y, z] = size(tnsr);
    t2 = reshape(tnsr, [x, y*z]);
    t3 = zscore(t2); % z-transformed
    %% SVD Analysis
    [~,s,~] = svds(t3,x); % We don't need bigger than x size
    sv = diag(s);
    pdet = prod(sv(sv ~= 0)); % pseudo determinant due to matrix size
    score = 0.5 * ((y*z) * log(2 * pi * exp(1)) + log( pdet));
end