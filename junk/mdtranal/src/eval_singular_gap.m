function [score] = eval_singular_gap(tnsr)
    [x, y, z] = size(tnsr);
    t2 = reshape(tnsr, [x, y*z]);
    %% SVD Analysis
    [~,s,~] = svds(t2,x); % We don't need bigger than x size
    sv = diag(s);
    %sv = sv / sum(abs(sv));
    score = sum(sv(1:x-1) - sv(2:x));
end