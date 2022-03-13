function [p, q, conv_flag] = power_fixed_point_ul(H, gamma, eta, p_ini, n_iter, quant_mode, decode_mode)

p_before = p_ini;

diff_norm = zeros(n_iter, 1);

for i = 1:n_iter
    
    if quant_mode == false
        q_before = q_from_p_ul_ind(H, p_before, eta);
    else
        q_before = q_from_p_ul_wz(H, p_before, eta);
    end
    
    p_after = p_update_ul(H, p_before, q_before, gamma, decode_mode);

    diff_norm(i) = norm(p_after - p_before);
    
    p_before = p_after;
end

p = p_after;

if quant_mode == false
    q = q_from_p_ul_ind(H, p_after, eta);
else
    q = q_from_p_ul_wz(H, p_after, eta);
end

conv_flag = false;
if (diff_norm(n_iter) - diff_norm(n_iter-1)) > 1
    conv_flag = true;
end

plot(1:n_iter, 10 * log10(diff_norm))
grid on
xlabel('iteartion number (n)')
ylabel('||p_{n+1} - p_n|| (dB)')

end