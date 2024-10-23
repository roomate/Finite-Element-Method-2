function grad = grad_g(r)

grad(1) = 2/3*r(1)^(-1/3)*sin(2*r(2)/3);
grad(2) = 2/3*r(1)^(-1/3)*cos(2*r(2)/3);