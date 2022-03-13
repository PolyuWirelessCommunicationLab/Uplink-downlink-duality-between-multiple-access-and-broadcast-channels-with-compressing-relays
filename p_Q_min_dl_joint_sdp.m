function [min_fin] = p_Q_min_dl_joint_sdp(H, gamma, eta)

[M, K] = size(H);

I_M = eye(M);
   
cvx_begin
    cvx_precision high
    
    for k = 1:K
        variable V
    end
    variable p(K)
    variable Q(M, M) symmetric
    minimize ( sum(p) + trace(Q) )
    subject to
    
    for k = 1:K
       H(:, k)' * V * diag(p) * V' * H(:, k) ...
           + H(:, k)' * Q * H(:, k) <= ...
                (1 + 1/gamma(k)) * ...
                    H(:, k)' * V(:, k) * p(k) * V(:, k)' * H(:, k);
    end
    
    eta(M)/(eta(M) + 1) * Q(M, M) ... 
        - 1/(eta(M) + 1) * I_M(:, M)' * V * diag(p) * V' * I_M(:, M) >= 0;
    
    for m = (M-1):-1:1   
         [eta(m)/(eta(m) + 1) * Q(m, m) ...
             - 1/(eta(m) + 1) * ...
                I_M(:, m)' * V * diag(p) * V' * I_M(:, m), ...
                    Q(m, (m+1):M); ...
          Q((m+1):M, m), Q((m+1):M, (m+1):M)] >= 0;
    end
          
    p >= zeros(K,1);
    Q >= 0;
cvx_end
 
p_fin = p;
Q_fin = Q;
min_fin = cvx_optval;
end