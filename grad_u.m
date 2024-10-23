function g_uf = grad_u(x)
g_uf(1) = pi*cos(pi*x(1))*sin(pi*x(2));
g_uf(2) = pi*sin(pi*x(1))*cos(pi*x(2));
