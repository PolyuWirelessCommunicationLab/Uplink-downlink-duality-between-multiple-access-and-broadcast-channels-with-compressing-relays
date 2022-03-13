function [p1_fin, p2_fin] = power_fixed_point_ul_m2_k2(p1_ini, p2_ini)

p1_ul = p1_ini;
p2_ul = p2_ini;

for i = 1:1000
    p1_ul = gamma1 * 1/(h1' * ((p2_ul * h2 * h2' + eye(3)) \ h1));
    p2_ul = gamma2 * 1/(h2' * ((p1_ul * h1 * h1' + eye(3)) \ h2));
end

p1_fin = p1_ul;
p2_fin = p2_ul;

end

