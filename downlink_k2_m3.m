H = normrnd(0, 1, 3, 2);

h1 = H(:,1);
h2 = H(:,2);

V = normrnd(0, 1, 3, 2);

v1_ini = V(:,1)/norm(V(:,1));
v2_ini = V(:,2)/norm(V(:,2));

gamma1 = (v1_ini' * h1)^2 / ((v2_ini' * h1)^2 + 1);
gamma2 = (v2_ini' * h2)^2 / ((v1_ini' * h2)^2 + 1);

A1 = [h1' zeros(1, 3); zeros(1, 3) h1'];
b1 = [sqrt(1 + 1/gamma1) * h1' zeros(1, 3)];
A2 = [h2' zeros(1, 3); zeros(1, 3) h2'];
b2 = [zeros(1, 3) sqrt(1 + 1/gamma2) * h2'];

cvx_begin
    variable v(6)
    minimize ( square_pos(norm(v(1:3))) + square_pos(norm(v(4:6))) )
    subject to
        norm([A1 * v; 1]) <= b1 * v
        norm([A2 * v; 1]) <= b2 * v
cvx_end

cvx_begin
    variable v(6)
    variable sqrtp(2)
    minimize ( sqrtp(1) + sqrtp(2) )
    subject to
        norm([A1 * v; 1]) <= b1 * v
        norm([A2 * v; 1]) <= b2 * v
        norm(v(1:3)) <= sqrtp(1)
        norm(v(4:6)) <= sqrtp(2)
cvx_end

v1_sol = v(1:3);
v2_sol = v(4:6);

gamma1_sol = (v1_sol' * h1)^2 / ((v2_sol' * h1)^2 + 1);
gamma2_sol = (v2_sol' * h2)^2 / ((v1_sol' * h2)^2 + 1);

tw1 = (p2_fin * h2 * h2' + eye(3)) \ h1;
tw2 = (p1_fin * h1 * h1' + eye(3)) \ h2;

w1 = tw1/norm(tw1);
w2 = tw2/norm(tw2);

cvx_begin
    variable p_dl(2)
    minimize (p_dl(1) + p_dl(2))
    subject to
        p_dl(2) * (w2' * h1)^2 + 1 <= 1/gamma1 * p_dl(1) * (w1' * h1)^2
        p_dl(1) * (w1' * h2)^2 + 1 <= 1/gamma2 * p_dl(2) * (w2' * h2)^2
cvx_end

p1_dl = p_dl(1);
p2_dl = p_dl(2);

w1_sol = w1;
w2_sol = w2;

gamma1_sol_ul = p1_dl * (w1_sol' * h1)^2 / (p2_dl * (w2_sol' * h1)^2 + 1);
gamma2_sol_ul = p2_dl * (w2_sol' * h2)^2 / (p1_dl * (w1_sol' * h2)^2 + 1);