function  [L,erreur_energie] = grad_u_uh(Numtri,Coorneu,XXF,face,Nbtri, UUh)
erreur_energie = 0;
L = zeros(Nbtri,1);
for l=1:Nbtri %On calcule le gradient brisé sur chacun des simplex
    grad_uh = 0;
    sommet1 = Numtri(l,1);
    sommet2 = Numtri(l,2);
    sommet3 = Numtri(l,3);
    
    Face1 = [sommet1 sommet2]; %Numéro des sommets de la face 1
    Face2 = [sommet2 sommet3]; %Numéro des sommets de la face 2
    Face3 = [sommet1 sommet3]; %Numéro des sommets de la face 3
    
    AB = [Coorneu(sommet2,1) - Coorneu(sommet1,1), Coorneu(sommet2,2) - Coorneu(sommet1,2)];
    AC = [Coorneu(sommet3,1) - Coorneu(sommet1,1), Coorneu(sommet3,2) - Coorneu(sommet1,2)];
    K = abs(AB(1)*AC(2) - AB(2)*AC(1))/2;
    
%     %On s'occupe de la face 1
%    
    i = find_index(Face1,face);  
    XXF1 = XXF(i,:);   
    UUh1 = UUh(i);
%     F1 = norm(Coorneu(sommet1,:) - Coorneu(sommet2,:)); 
%     X = [XXF1 - Coorneu(sommet3,:)];
%     n1 = [ - (Coorneu(sommet2,2) - Coorneu(sommet1,2)), Coorneu(sommet2,1) - Coorneu(sommet1,1)];
%     if (scalaire(n1,X)<0)
%         n1 = -n1;
%     end
%      n1 = n1/norm(n1);
%      grad_uh = F1*UUh1*n1/K; %La somme porte sur les 3 faces de K
%     
%     
%     %On s'occupe de la face 2
%    
    j = find_index(Face2,face);
    XXF2 = XXF(j,:);
    UUh2 = UUh(j);
%     X = [XXF2 - Coorneu(sommet1,:)];
%     n2 = [ - (Coorneu(sommet3,2) - Coorneu(sommet2,2)), Coorneu(sommet3,1) - Coorneu(sommet2,1)];
%     F2 = norm(Coorneu(sommet3,:) - Coorneu(sommet2,:));
%     if (scalaire(n2,X)<0)
%          n2 = - n2;
%     end
%     n2 = n2/norm(n2);
%     grad_uh = grad_uh + F2*UUh2/K*n2;
%     
%     
%     %On s'occupe de la face 3
%         
    p = find_index(Face3, face);
    XXF3 = XXF(p,:);
%     X = [XXF3 - Coorneu(sommet2,:)];
    UUh3 = UUh(p);
%     F3 = norm(Coorneu(sommet3,:) - Coorneu(sommet1,:));
%     n3 = [ - (Coorneu(sommet3,2) - Coorneu(sommet1,2)), Coorneu(sommet3,1) - Coorneu(sommet1,1)];
%     if (scalaire(n3,X)<0)
%        n3 = - n3;
%     end
%     n3 = n3/norm(n3);
%     grad_uh = grad_uh + F3*UUh3/K*n3;

    A = [XXF1(1), XXF1(2), 1; XXF2(1), XXF2(2), 1; XXF3(1), XXF3(2), 1];
    U = [UUh1; UUh2; UUh3];
    Coeff = A\U;
    
    
    
    %%
    grad_u_uh1 = norm(grad_u(XXF2) - Coeff(1:2))^2 +  norm(grad_u(XXF3) - Coeff(1:2))^2 + norm(grad_u(XXF1) - Coeff(1:2))^2; 
    
    L(l) = grad_u_uh1*K/3;
    erreur_energie = erreur_energie + K/3*grad_u_uh1; %On somme sur l'ensemble des mailles
end
erreur_energie = sqrt(erreur_energie);
    