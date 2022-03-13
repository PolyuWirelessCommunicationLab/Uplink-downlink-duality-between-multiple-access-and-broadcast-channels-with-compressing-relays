function [q] = q_from_p_ul_wz(H, p, eta)

[M, ~] = size(H);

q = zeros(M, 1);

q(1) = 1/eta(1) * (H(1, :) * diag(p) * H(1, :)' + 1);

for m = 2:M
    Sigma = H * diag(p) * H' + diag(q) + eye(M);
    q(m) = 1/eta(m) * (H(m, :) * diag(p) * H(m, :)' + 1 ...
        - Sigma(m,1:(m-1)) * (Sigma(1:(m-1), 1:(m-1)) \ Sigma(1:(m-1),m)));
end