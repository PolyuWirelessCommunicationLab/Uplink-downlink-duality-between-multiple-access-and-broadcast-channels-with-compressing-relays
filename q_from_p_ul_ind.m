function [q] = q_from_p_ul_ind(H, p, eta)

[M, ~] = size(H);

q = zeros(M,1);

for m = 1:M
    q(m) = 1/eta(m) * (H(m, :) * diag(p) * H(m, :)' + 1);
end