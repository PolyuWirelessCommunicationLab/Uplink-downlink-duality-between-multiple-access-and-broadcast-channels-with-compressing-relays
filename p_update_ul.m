function [p] = p_update_ul(H, p, q, gamma, decode_mode)

[M, K] = size(H);

if decode_mode == false
    for k = 1:K
        inter_mat = ...
            H * diag(p) * H' - H(:,k) * p(k) * H(:,k)' + diag(q) + eye(M);
        p(k) = gamma(k) * 1/(H(:, k)' * (inter_mat \ H(:, k)));
    end
else
    for k = 1:(K-1)
        inter_mat = ...
           H(:, (k+1):K) * diag(p((k+1):K)) * H(:, (k+1):K)'  + diag(q) + eye(M);
        p(k) = gamma(k) * 1/(H(:, k)' * (inter_mat \ H(:, k)));
    end
    inter_mat = ...
        0 + diag(q) + eye(M);
    p(K) = gamma(K) * 1/(H(:, K)' * (inter_mat \ H(:, K)));    
end
    
end