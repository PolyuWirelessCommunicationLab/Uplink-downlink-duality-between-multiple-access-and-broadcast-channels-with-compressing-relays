M = 3;
K = 3;

seed = 189;
rng(seed)

H = normrnd(0, 1, M, K);

p_ini = 0.1 * ones(K,1);

n_iter = 100;

r_range = 1.5:0.05:3;

p_ul_ind_tin_range = zeros(length(r_range), 1);

p_ul_ind_sic_range = zeros(length(r_range), 1);

p_ul_wz_tin_range = zeros(length(r_range), 1);

p_ul_wz_sic_range = zeros(length(r_range), 1);

p_dl_ind_tin_range = zeros(length(r_range), 1);

p_dl_ind_dpc_range = zeros(length(r_range), 1);

p_dl_mc_tin_range = zeros(length(r_range), 1);

p_dl_mc_dpc_range = zeros(length(r_range), 1);

c = 3 * ones(M, 1);

eta = 2.^(2 * c) - 1;

for i = 1:length(r_range)
   
    r = r_range(i) * ones(K, 1);
    
    gamma = 2.^(2 * r) - 1;
    
    [p, q, conv_flag] = power_fixed_point_ul(H, gamma, eta, p_ini, n_iter, false, false);
    
    p_ul_ind_tin_range(i) = sum(p);
    
    W = mmse_beamformer_ul(H, p, q, false);

    [~, ~, min_fin] = p_Q_min_dl(H, gamma, eta, W, false, false);
    
    p_dl_ind_tin_range(i) = min_fin;
    
    if conv_flag == true
        p_ul_ind_tin_range(i) = Inf;
    end
    
    [p, q, conv_flag] = power_fixed_point_ul(H, gamma, eta, p_ini, n_iter, false, true);
    
    p_ul_ind_sic_range(i) = sum(p);
    
    W = mmse_beamformer_ul(H, p, q, true);

    [~, ~, min_fin] = p_Q_min_dl(H, gamma, eta, W, false, true);
    
    p_dl_ind_dpc_range(i) = min_fin;
    
    if conv_flag == true
        p_ul_ind_sic_range(i) = Inf;
    end

    [p, q, conv_flag] = power_fixed_point_ul(H, gamma, eta, p_ini, n_iter, true, false);
    
    p_ul_wz_tin_range(i) = sum(p);
    
    W = mmse_beamformer_ul(H, p, q, false);
    
    [~, ~, min_fin] = p_Q_min_dl(H, gamma, eta, W, true, false);
    
    p_dl_mc_tin_range(i) = min_fin;
    
    if conv_flag == true
        p_ul_wz_tin_range(i) = Inf;
    end
    
    [p, q, conv_flag] = power_fixed_point_ul(H, gamma, eta, p_ini, n_iter, true, true);
    
    p_ul_wz_sic_range(i) = sum(p);
    
    W = mmse_beamformer_ul(H, p, q, true);
    
    [~, ~, min_fin] = p_Q_min_dl(H, gamma, eta, W, true, true);
    
    p_dl_mc_dpc_range(i) = min_fin;
    
    if conv_flag == true
        p_ul_wz_sic_range(i) = Inf;
    end
end

plot(r_range, 10 * log10(p_ul_ind_tin_range), '-')
hold on
grid on
plot(r_range, 10 * log10(p_ul_ind_sic_range), '--')
plot(r_range, 10 * log10(p_ul_wz_tin_range), ':')
plot(r_range, 10 * log10(p_ul_wz_sic_range), '-.')
plot(r_range, 10 * log10(p_dl_ind_tin_range), 'o')
plot(r_range, 10 * log10(p_dl_ind_dpc_range), '+')
plot(r_range, 10 * log10(p_dl_mc_tin_range), 'x')
plot(r_range, 10 * log10(p_dl_mc_dpc_range), '*')
legend('UL (IN, TIN)', 'UL (IN, SIC)', 'UL (WZ, TIN)', 'UL (WZ, SIC)', ...
    'DL (IN, LIN)', 'DL (IN, DPC)', 'DL (MV, LIN)', 'DL (MV, DPC)', ...
    'Location', 'best')
xlabel('Individual user rate target (bps)')
ylabel('Sum power (dB)')
xlim([1.4 3])