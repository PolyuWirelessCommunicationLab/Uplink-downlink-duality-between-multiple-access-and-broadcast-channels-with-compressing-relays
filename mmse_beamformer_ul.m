function W = mmse_beamformer_ul(H, p, q, decode_mode)

[M, K] = size(H);

W = zeros(M, K);

if decode_mode == false
    for k = 1:K
        inter_mat = H * diag(p) * H' - H(:,k) * p(k) * H(:,k)' + diag(q) + eye(M);
        w = inter_mat \ H(:, k);
        W(:,k) = w / norm(w);
    end
else
    for k = 1:(K-1)
        inter_mat = H(:, (k+1):K) * diag(p((k+1):K)) * H(:, (k+1):K)' + diag(q) + eye(M);
        w = inter_mat \ H(:, k);
        W(:,k) = w / norm(w);
    end    
        inter_mat = 0 + diag(q) + eye(M);
        w = inter_mat \ H(:, K);
        W(:,K) = w / norm(w);    
end

end