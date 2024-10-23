function affiche2(UU,Numtri,Coorneu,Numface,XXF, titre)
Nbpt = size(Coorneu,1);
Nbtri = size(Numtri,1);
for i = 1:Nbtri
    sommet1 = Numtri(i,1);
    sommet2 = Numtri(i,2);
    sommet3 = Numtri(i,3);
    
    m = Numface(i,1);
    j = Numface(i,2);
    p = Numface(i,3);
   
    
    U = [UU(m); UU(j); UU(p)];
    A = [XXF(m,1), XXF(m,2), 1; XXF(j,1), XXF(j,2), 1; XXF(p,1), XXF(p,2), 1];
    Coeff = A\U;

   Usommet1 = scalaire(Coeff,[Coorneu(sommet1,1); Coorneu(sommet1,2); 1]); 
   Usommet2 = scalaire(Coeff,[Coorneu(sommet2,1); Coorneu(sommet2,2); 1]);
   Usommet3 = scalaire(Coeff,[Coorneu(sommet3,1); Coorneu(sommet3,2); 1]); 

   U(1) =  Usommet1;
   U(2) =  Usommet2;
   U(3) =  Usommet3;
   
   X = [Coorneu(sommet1,1) Coorneu(sommet2,1) Coorneu(sommet3,1)];
   Y = [Coorneu(sommet1,2) Coorneu(sommet2,2) Coorneu(sommet3,2)];
   Z = [U(1), U(2), U(3)];

   patch(X,Y,Z, [0, 3/5, 0]);
   view(3);
end
title(titre);
    
    
    