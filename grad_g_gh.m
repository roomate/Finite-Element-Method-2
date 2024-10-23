function [L,erreur_energie] = grad_g_gh(Numtri,Coorneu,XXF,face,Nbtri, Ugh)
erreur_energie = 0;
L = zeros(Nbtri,1);
for l=1:Nbtri %On calcule le gradient brisé sur chacun des simplex
    grad_ugh = 0;
    grad_u = 0;
    sommet1 = Numtri(l,1);
    sommet2 = Numtri(l,2);
    sommet3 = Numtri(l,3);
    
    Face1 = [sommet1 sommet2]; %Numéro des sommets de la face 1
    Face2 = [sommet2 sommet3]; %Numéro des sommets de la face 2
    Face3 = [sommet1 sommet3]; %Numéro des sommets de la face 3
    
    AB = [Coorneu(sommet2,1) - Coorneu(sommet1,1), Coorneu(sommet2,2) - Coorneu(sommet1,2)];
    AC = [Coorneu(sommet3,1) - Coorneu(sommet1,1), Coorneu(sommet3,2) - Coorneu(sommet1,2)];
    K = abs(AB(1)*AC(2) - AB(2)*AC(1))/2;
    
    i = find_index(Face1,face);  
    XXF1 = XXF(i,:);   
    Ugh1 = Ugh(i);
    
    j = find_index(Face2,face);
    XXF2 = XXF(j,:);
    Ugh2 = Ugh(j);
    
    p = find_index(Face3, face);
    XXF3 = XXF(p,:);
    Ugh3 = Ugh(p);
    
    %% On s'occupe de la face 1
   
    F1 = norm(Coorneu(sommet1,:) - Coorneu(sommet2,:)); 
    X = [XXF1 - Coorneu(sommet3,:)];
    n1 = [ - (Coorneu(sommet2,2) - Coorneu(sommet1,2)), Coorneu(sommet2,1) - Coorneu(sommet1,1)];
    if (scalaire(n1,X)<0)
        n1 = -n1;
    end
     n1 = n1/norm(n1);
     grad_ugh = F1*Ugh1*n1/K; %La somme porte sur les 3 faces de K
     grad_u = F1/K*g(polar(XXF1))*n1;
    
    %% On s'occupe de la face 2
   
    X = [XXF2 - Coorneu(sommet1,:)];
    n2 = [ - (Coorneu(sommet3,2) - Coorneu(sommet2,2)), Coorneu(sommet3,1) - Coorneu(sommet2,1)];
    F2 = norm(Coorneu(sommet3,:) - Coorneu(sommet2,:));
    if (scalaire(n2,X)<0)
         n2 = -n2;
    end
    n2 = n2/norm(n2);
    grad_ugh = grad_ugh + F2*Ugh2/K*n2;
    grad_u = grad_u + F2/K*g(polar(XXF2))*n2;
    
    %% On s'occupe de la face 3
    X = [XXF3 - Coorneu(sommet2,:)];
    F3 = norm(Coorneu(sommet3,:) - Coorneu(sommet1,:));
    n3 = [ - (Coorneu(sommet3,2) - Coorneu(sommet1,2)), Coorneu(sommet3,1) - Coorneu(sommet1,1)];
    if (scalaire(n3,X)<0)
       n3 = - n3;
    end
    n3 = n3/norm(n3);
    grad_ugh = grad_ugh + F3*Ugh3/K*n3;
    grad_u = grad_u + F3/K*n3*g(polar(XXF3));

    %%
    polar1 = polar(XXF1);
    polar2 = polar(XXF2);
    polar3 = polar(XXF3);
    grad_u_uh1 = norm(cart(grad_g(polar2),polar2(2)) - grad_ugh)^2 +  norm(cart(grad_g(polar3),polar3(2)) - grad_ugh)^2 + norm(cart(grad_g(polar1),polar1(2)) - grad_ugh)^2;
    L(l) = K/3*grad_u_uh1;
    erreur_energie = erreur_energie + K/3*grad_u_uh1; %On somme sur l'ensemble des mailles
end
erreur_energie = sqrt(erreur_energie);