function [p_fin, Q_fin, min_fin] = p_Q_min_dl(H, gamma, eta, W, quant_mode, encode_mode)

[M, K] = size(H);

V = W;

I_M = eye(M);
   
cvx_begin sdp
    cvx_solver sedumi
    cvx_precision high
    
    variable p(K)
    variable term(K)
    
    if quant_mode == true
        variable Q(M, M) symmetric
    else
        variable Q(M, M) diagonal
    end
    minimize ( sum(p) + trace(Q) )
    subject to
    
    if encode_mode == false        
        for k = 1:K
           H(:, k)' * V * diag(p) * V' * H(:, k) ...
               + H(:, k)' * Q * H(:, k) + 1 <= ...
                    (1 + 1/gamma(k)) * ...
                        H(:, k)' * V(:, k) * p(k) * V(:, k)' * H(:, k);
        end
    else
        0 + H(:, 1)' * Q * H(:, 1) + 1 <= ...
                    1/gamma(1) * ...
                        H(:, 1)' * V(:, 1) * p(1) * V(:, 1)' * H(:, 1);         
        for k = 2:K
            H(:, k)' * V(:, 1:(k-1)) * diag(p(1:(k-1))) * V(:, 1:(k-1))' * H(:, k) ...
               + H(:, k)' * Q * H(:, k) + 1 <= ...
                    1/gamma(k) * ...
                        H(:, k)' * V(:, k) * p(k) * V(:, k)' * H(:, k);
        end
    end
    
    Q(M, M) - 1/eta(M) * I_M(:, M)' * V * diag(p) * V' * I_M(:, M) >= 0;
    
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