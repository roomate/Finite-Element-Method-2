function [a,b] = reg_lin(X,Y)

Z = [ones(length(X),1), X];
w = Z\Y;
b = w(1);
a = w(2);