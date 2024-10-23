function produit_scalaire = scalaire(A,B)
produit_scalaire = 0;
if (size(A) ~= size(B))
    disp("Les vecteurs ne sont pas de même dimensions");
end
for i=1:max(size(A,1),size(A,2))
    produit_scalaire = produit_scalaire + A(i)*B(i);
end