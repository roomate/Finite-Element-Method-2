% =====================================================
% principal_stationnaire;
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour :
%
% 1) l'equation suivante stationnaire, avec conditions de
% Dirichlet homogene
%  -laplacien(u) = f
%
%
% =====================================================

% Donnees du probleme.
nom_maillage = 'geomRect.msh';
affichage = true; % false; %

% Lecture du maillage et affichage.
[Nbpt, Nbtri, Coorneu, Refneu, Numtri, Reftri,Nbaretes,Numaretes,Refaretes] = lecture_msh2(nom_maillage);
%Comptage du nombre de face intérieur pour avoir la taille des matrices
Nbface_int = 3/2*Nbtri - Nbaretes/2;

% Declarations des matrices EF et vecteur second membre.
KK = sparse(Nbface_int,Nbface_int); %Matrice de rigidité
B = zeros(Nbface_int,1); %Second membre pour 1er cas
Bg = zeros(Nbface_int, Nbaretes); % Second membre pour le second cas
G = zeros(Nbaretes,1);
XF = zeros(Nbface_int,2); %Liste des coordonnées du milieu de chacune des faces
GG = zeros(Nbpt,1);
UU = zeros(Nbpt,1);

for i=1:Nbpt
    x = Coorneu(i,1);
    y = Coorneu(i,2);
    if (x == 0) 
        GG(i) = 0;
    else
        GG(i) = g(polar([x,y]));
    end
end

for i=1:Nbpt
    x = Coorneu(i,1);
    y = Coorneu(i,2);
    UU(i) = f([x,y])/(2*pi^2);
end

%Faire une liste des face intérieurs
face_int = zeros(Nbface_int,2);
Numface = zeros(Nbtri,3); %Numéro des 3 faces dans la maille
face = zeros(Nbface_int + Nbaretes,2);
k = 1;
for n=1:Nbtri
    sommet1 = Numtri(n,1);
    sommet2 = Numtri(n,2);
    sommet3 = Numtri(n,3);
    if ((membre([sommet1,sommet2],face_int) == false) && (membre([sommet1, sommet2],Numaretes) == false) )
        face_int(k,:) = [sommet1, sommet2];
        face(k,:) = face_int(k,:);        
        k = k + 1;
    end
    if ((membre([sommet1,sommet3],face_int) == false) && (membre([sommet1, sommet3],Numaretes) == false))
        face_int(k,:) = [sommet1, sommet3];
        face(k,:) = face_int(k,:);
        k = k + 1;
    end
    if ((membre([sommet2,sommet3],face_int) == false) && (membre([sommet3, sommet2],Numaretes) == false))
        face_int(k,:) = [sommet2, sommet3];
        face(k,:) = face_int(k,:);
        k = k + 1;
    end
end

%On complète la liste face
for i=1:Nbaretes
    face(Nbface_int + i,:) = Numaretes(i,:);

end
%On a donc une 1re partie avec face intérieur puis que les aretes




%%
% Boucle d'assemblage de K et B
for l=1:Nbtri
    sommet1 = Numtri(l,1);
    sommet2 = Numtri(l,2);
    sommet3 = Numtri(l,3);
    
    Face1 = [sommet1 sommet2]; %Numéro des sommets de la face 1
    Face2 = [sommet2 sommet3]; %Numéro des sommets de la face 2
    Face3 = [sommet1 sommet3]; %Numéro des sommets de la face 3
    
    %On met toutes les faces dans Numface, pour avoir toutes les faces de
    %tous les triangles
    Numface(l,1) = find_index(Face1,face); %Numéro de la Face1
    Numface(l,2) = find_index(Face2,face); %Numéro de la Face2
    Numface(l,3) = find_index(Face3,face); %Numéro de la Face3
    
    %Mesure du triangle
    AB = [Coorneu(sommet2,1) - Coorneu(sommet1,1), Coorneu(sommet2,2) - Coorneu(sommet1,2)]; %Face 1
    AC = [Coorneu(sommet3,1) - Coorneu(sommet1,1), Coorneu(sommet3,2) - Coorneu(sommet1,2)]; % Face 3
    BC = [Coorneu(sommet3,1) - Coorneu(sommet2,1), Coorneu(sommet3,2) - Coorneu(sommet2,2)]; %Face 2
    K = abs(AB(1)*AC(2) - AB(2)*AC(1))/2;
    
    %Calcul des normales des faces
    
    xF1 = [ Coorneu(sommet1,1) + Coorneu(sommet2,1) , Coorneu(sommet1,2) + Coorneu(sommet2,2)]/2; %Milieu de la face1
    X = [xF1 - Coorneu(sommet3,:)];
    n1 = [ - (Coorneu(sommet2,2) - Coorneu(sommet1,2)), Coorneu(sommet2,1) - Coorneu(sommet1,1)];
    n1 = n1/norm(n1);
    if (scalaire(n1,X)<0)
        n1 = -n1;
    end
    
    xF2 = [ Coorneu(sommet3,1) + Coorneu(sommet2,1) , Coorneu(sommet3,2) + Coorneu(sommet2,2)]/2; %Milieu de la face2
    X = [xF2 - Coorneu(sommet1,:)];
    n2 = [ - (Coorneu(sommet3,2) - Coorneu(sommet2,2)), Coorneu(sommet3,1) - Coorneu(sommet2,1)];
    n2 = n2/norm(n2);
    F2 = norm(Coorneu(sommet3,:) - Coorneu(sommet2,:));
    if (scalaire(n2,X)<0)
        n2 = -n2;
    end
    
    xF3 = [ Coorneu(sommet3,1) + Coorneu(sommet1,1) , Coorneu(sommet1,2) + Coorneu(sommet3,2)]/2; %Milieu de la face3
    X = [xF3 - Coorneu(sommet2,:)];
    F3 = norm(Coorneu(sommet3,:) - Coorneu(sommet1,:));
    n3 = [ - (Coorneu(sommet3,2) - Coorneu(sommet1,2)), Coorneu(sommet3,1) - Coorneu(sommet1,1)];
    if (scalaire(n3,X)<0)
        n3 = - n3;
    end
    n3 = n3/norm(n3);
    %Index des faces respectives dans la liste face SI la face n'est pas une arête
    bool1 = membre(Face1,face_int); %Vrai si la face 1 n'est pas une arête
    if (bool1 == true)
            i = find_index(Face1,face_int);
            XF(i,:) = xF1;
            F1 = norm(Coorneu(sommet1,:) - Coorneu(sommet2,:)); 
            KK(i,i) = KK(i,i) + F1^2/K; %Il y a deux termes à sommer sur la diagonale
            B(i) = B(i) + f(xF1)/3*K;    
    end
    bool2 = membre(Face2,face);
    if (bool2 == true)
            j = find_index(Face2,face_int);
            XF(j,:) = xF2;
            KK(j,j) = KK(j,j) + F2^2/K; %Il y a deux termes à sommer sur la diagonale
            B(j) = B(j) + f(xF2)/3*K;
    end
    bool3 = membre(Face3,face_int);
    if (bool3 == true)
        p = find_index(Face3,face_int);
        XF(p,:) = xF3;
        KK(p,p) = KK(p,p) + F3^2/K; %Il y a deux termes à sommer sur la diagonale
        B(p) = B(p) + f(xF3)/3*K;
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%Assemblage de la matrice de rigidité KK et du second membre B%%%%%%%%
    if (bool1 == true && bool2 == true)
        KK(i,j) = F1*F2/K*scalaire(n1,n2); %Il n'y en a qu'un seul sur les termes non diagonaux
        KK(j,i) = KK(i,j);
    end
    if (bool1 == true && bool3 == true)
        KK(i,p) = F1*F3/K*scalaire(n1,n3);
        KK(p,i) = KK(i,p);
    end
    if (bool2 == true && bool3 == true)
        KK(j,p) = F2*F3/K*scalaire(n2,n3);
        KK(p,j) =  KK(j,p); 
    end
end
%%
XXF = zeros(Nbface_int + Nbaretes,2);
XXF(1:Nbface_int,1) = XF(:,1);
XXF(1:Nbface_int,2) = XF(:,2);
for i = 1:Nbaretes
    XXF(Nbface_int + i,:) = [Coorneu(Numaretes(i,1),1) + Coorneu(Numaretes(i,2),1), Coorneu(Numaretes(i,1),2) + Coorneu(Numaretes(i,2),2)]/2; 
end
for i = 1:Nbaretes
    G(i) = g(polar(XXF(Nbface_int + i,:)));
end
%%
%Boucle d'assemblage de Bg
for i=1:Nbaretes
    sommet1 = Numaretes(i,1);
    sommet2 = Numaretes(i,2);
    
    l = find_index2(i + Nbface_int,Numface); %Donne le numéro de triangle comportant cette face extérieur
    k = 1;
    while (Numface(l,k) ~= i + Nbface_int)
        k = k + 1;
    end
    NumFace2 = Numface(l,mod(k+1,3) + 1);
    NumFace3 = Numface(l,mod(k+2,3) + 1);
    Face2 = face(NumFace2,:);
    
    if (Face2(1) == sommet1)
        sommet3 = Face2(2);
    else if (Face2(1) == sommet2)
        sommet3 = Face2(2);
        else
            sommet3 = Face2(1);
        end
    end
    Face_exte = [sommet1, sommet2];
    Face2 = [sommet1, sommet3];
    Face3 = [sommet2, sommet3];
    F1 = norm(Coorneu(sommet1,:) - Coorneu(sommet2,:));
    j = find_index(Face2,face);
    p = find_index(Face3, face);
    Face_exte = [Coorneu(sommet1,:) - Coorneu(sommet2,:)];
    BC = [Coorneu(Face2(1),:) - Coorneu(Face2(2),:)];
    AC = [Coorneu(Face3(1),:) - Coorneu(Face3(2),:)];
    F2 = norm(Coorneu(sommet3,:) - Coorneu(sommet1,:));
    F3 = norm(Coorneu(sommet3,:) - Coorneu(sommet2,:));
    K = abs(AC(1)*BC(2) - AC(2)*BC(1))/2;
    xF1 = XXF(Nbface_int + i,:);
    xF2 = XXF(NumFace2,:);
    xF3 = XXF(NumFace3,:);
    
    X = [xF1 - Coorneu(sommet3,:)];
    n1 = [ - (Coorneu(sommet2,2) - Coorneu(sommet1,2)), Coorneu(sommet2,1) - Coorneu(sommet1,1)];
    n1 = n1/norm(n1);
    if (scalaire(n1,X)<0)
        n1 = -n1;
    end
    
    X = [xF2 - Coorneu(sommet2,:)];
    n2 = [ - (Coorneu(sommet3,2) - Coorneu(sommet1,2)), Coorneu(sommet3,1) - Coorneu(sommet1,1)];
    n2 = n2/norm(n2);
    if (scalaire(n2,X)<0)
        n2 = -n2;
    end
    
    X = [xF3 - Coorneu(sommet1,:)];
    n3 = [ - (Coorneu(sommet3,2) - Coorneu(sommet2,2)), Coorneu(sommet3,1) - Coorneu(sommet2,1)];
    n3 = n3/norm(n3);
    if (scalaire(n3,X)<0)
        n3 = -n3;
    end
    
    Bg(j,i) =  Bg(j,i) + F1*F2/K*scalaire(n1,n2);
    Bg(p,i) = Bg(p,i) + F1*F3/K*scalaire(n1,n3);
end
    
    
%%
% Résolution du problème par inversion et sortie des résultats
%validation et résolution numérique.
Uh = KK\B;

tUh = -(KK\Bg)*G;

UUh = zeros(Nbface_int + Nbaretes,1);
Ugh = zeros(Nbface_int + Nbaretes,1);


UUh(1:Nbface_int) = Uh; %On met à 0 les termes sur les aretes
Ugh(1: Nbface_int) = tUh; %Pareil

for i = 1:Nbface_int + Nbaretes
    if (i <= Nbface_int)
    Ugh(i) = tUh(i);
    else
    Ugh(i) = g(polar(XXF(i,:)));
    end
end
%%
%Construction de l'erreur ||grad(u - uh)|| pour |Vh1|
[Lu,Grad_u_uh] = grad_u_uh(Numtri,Coorneu,XXF,face,Nbtri, UUh);  %Lu permettra de connaitre l'erreur dans chacune des mailles
GGrad_u_uh = [1.5477, 1.1790, 0.5835, 0.2910, 0.1454, 0.0727, 0.0363];
Vh1 = [4, 20, 88, 368, 1504, 6080, 24448];

%loglog(Vh1,GGrad_u_uh,"DisplayName", 'erreur en fonction de |Vh1|');
hold on;
GGrad_u_uh = GGrad_u_uh';
[a,b] = reg_lin(log(Vh1'), log(GGrad_u_uh));
Grad_u_uh2 = a*log(Vh1) + b;
%loglog(Vh1,exp(Grad_u_uh2), "DisplayName", "régression linéaire de l'erreur");
hold off
%%
%Construction de l'erreur ||grad(u - uh)|| pour h
GGrad_u_uh = [1.5477, 1.1790, 0.5835, 0.2910, 0.1454, 0.0727, 0.0363];
h = [1/2, 1/2^2, 1/2^3, 1/2^4, 1/2^5, 1/2^6, 1/2^7, 1/2^8];
%loglog(h,GGrad_u_uh);

%%
%Construction de l'erreur ||grad(g - gh)|| pour |Vh1|
[Lg,Grad_g_gh] = grad_g_gh(Numtri,Coorneu,XXF,face,Nbtri, Ugh); %Lg permettra de connaître l'erreur dans chacune des mailles
GGrad_g_gh = [0.3034,0.2094, 0.1387, 0.0898, 0.0575, 0.0366,0.0232];
Vh1 = [14, 64, 272, 1120, 4544, 18304, 73472];
%loglog(Vh1,GGrad_g_gh);
hold on;
[c, d] = reg_lin(log(Vh1'), log(GGrad_g_gh'));
Grad_g_gh2 = c*log(Vh1) + d;
%loglog(Vh1,exp(Grad_g_gh2));
%%
%Construction de l'erreur ||grad(g - gh)|| pour h
GGrad_g_gh = [0.3034,0.2094, 0.1387, 0.0898, 0.0575, 0.0366,0.0232];
h = [1/2, 1/2^2, 1/2^3, 1/2^4, 1/2^5, 1/2^6,1/2^7];
%loglog(h,GGrad_g_gh);
%%
%Reconstruction du potentiel (Exercice 2)
% %% 
% On utilise les vertex, qui sont les sommets des triangles pour construire
% ce potentiel.
Sh = zeros(Nbpt,1); %Initialisation de la liste Sh (cas 1)
Sgh = zeros(Nbpt,1); %Initialisation de la liste Sgh (cas 2)
iter = zeros(Nbpt,1); %Connaître le nombre de triangle connecté à un sommet intérieur
iter2 = zeros(Nbpt,1);
for i = 1:Nbface_int
    sommet1 = face(i,1); %On prend le 1er sommet de la face i
    sommet2 = face(i,2); %Et le second
    if (Refneu(sommet1) == 0) %Si le sommet n'est pas sur le bord
        Sh(sommet1) = Sh(sommet1) + Uh(i);
        iter(sommet1) = iter(sommet1) + 1;
    end
    if (Refneu(sommet2) == 0)
        Sh(sommet2) = Sh(sommet2) + Uh(i);
        iter(sommet2) = iter(sommet2) + 1;
    end
end

for i = 1:Nbface_int + Nbaretes
    sommet1 = face(i,1); %On prend le 1er sommet de la face i
    sommet2 = face(i,2); %Et le second
    Sgh(sommet1) = Sgh(sommet1) + Ugh(i);
    iter2(sommet1) = iter2(sommet1) + 1;
    Sgh(sommet2) = Sgh(sommet2) + Ugh(i);
    iter2(sommet2) = iter2(sommet2) + 1;
end


for i = 1:Nbpt
    if (iter(i) ~= 0)
        Sh(i) = Sh(i)/iter(i);
    end
    if (iter2(i) ~= 0)
        Sgh(i) = Sgh(i)/iter2(i);
    end
end

%%
%Exercice 4
%Estimation d'erreur


fh=zeros(Nbtri,1);
Xk=zeros(Nbtri,2);

for i=1:Nbtri
    S1=Numtri(i,1);       %sommets du triangle i
    S2=Numtri(i,2); 
    S3=Numtri(i,3);
    
    X1=[Coorneu(S1,1),Coorneu(S1,2)];   %Coordonnées des sommets du triangle i
    X2=[Coorneu(S2,1),Coorneu(S2,2)];
    X3=[Coorneu(S3,1),Coorneu(S3,2)];
    M1=(X1+X2)/2;                       %Milieux des côtés
    M2=(X1+X3)/2;
    M3=(X2+X3)/2;
    
    Y=f(M1)+f(M2)+f(M3);
    Y=Y/3;
    fh(i,1)=Y;
    
    xk=(X1+X2+X3)/3;                    %barycentre du triangle i  
    Xk(i,:)=[xk(1) xk(2)];
end

ETA_k2=zeros(Nbtri,1);

%Pour le cas1

for l=1:Nbtri
    
    i = Numface(l,1); 
    j = Numface(l,2);
    k = Numface(l,3);
    
    sommet1 = face(i,1);
    sommet2 = face(i,2);
    if (face(j,1) ~= sommet2 && face(j,1) ~=sommet1)
        sommet3 = face(j,1);
    else
        sommet3 = face(j,2);
    end
    Face1 = [sommet1, sommet2];
    Face2 = [sommet2, sommet3];
    Face3 = [sommet1, sommet3];
    
    j = find_index(Face2,face);
    k = find_index(Face3,face);
    
    XF1 = XXF(i,:);
    XF2 = XXF(j,:);
    XF3 = XXF(k,:);
    
    X1=[Coorneu(sommet1,1),Coorneu(sommet1,2)];   %Coordonnées des sommets du triangle i
    X2=[Coorneu(sommet2,1),Coorneu(sommet2,2)];
    X3=[Coorneu(sommet3,1),Coorneu(sommet3,2)];
    
    X = [XF1 - Coorneu(sommet3,:)];
    n1 = [ - (X2(2) - X1(2)), X2(1) - X1(1)];
    n1 = n1/norm(n1);
    if (scalaire(n1,X)<0)
        n1 = - n1;
    end
    
    X = [XF2 - Coorneu(sommet1,:)];
    n2 = [ - (X3(2) - X2(2)), X3(1) - X2(1)];
    n2 = n2/norm(n2);
    if (scalaire(n2,X)<0)
        n2 = - n2;
    end
    
    X = [XF3 - Coorneu(sommet2,:)];
    n3 = [ -( X3(2) - X1(2)), X3(1) - X1(1)];
    n3 = n3/norm(n3);
    if (scalaire(n3,X) < 0)
        n3 = - n3;
    end
    
    K=abs(det([X2-X1;X3-X1]))/2;             %Mesure de K
    F1 = norm(X2 - X1);
    F2 = norm(X3 - X2);
    F3 = norm(X3 - X1);
    
    xk=[Xk(l,1),Xk(l,2)];     %barycentre
    
    Sh1 = (Sh(sommet1) + Sh(sommet2))/2; %On prend Sh(xF) en moyennant sur les sommets adjacents de la face
    Sh2 = (Sh(sommet2) + Sh(sommet3))/2;
    Sh3 = (Sh(sommet3) + Sh(sommet1))/2;
    
    Uh1 = UUh(i);
    Uh2 = UUh(j);
    Uh3 = UUh(k);
    
    %Calcul du premier terme de l'estimateur d'erreur
    Y = (K/3)*fh(l,1)^2/4*(norm(XF1-xk)^2 + norm(XF2-xk)^2 + norm(XF3-xk)^2);
    
    %Calcul du second terme de l'estimateur d'erreur
    
    Y = Y + K * norm((Uh1 - Sh1)*F1/K*n1 + (Uh2 - Sh2)*F2/K*n2 + (Uh3 - Sh3)*F3/K*n3)^2;
    
    ETA_k2(l,1)=Y;   
end


%Pour le cas 2

ETA_k2_2 = zeros(Nbtri,1);

for l=1:Nbtri
    
    i = Numface(l,1); 
    j = Numface(l,2);
    k = Numface(l,3);
    
    sommet1 = face(i,1);
    sommet2 = face(i,2);
    
    if (face(j,1) ~= sommet2 && face(j,1) ~=sommet1)
        sommet3 = face(j,1);
    else
        sommet3 = face(j,2);
    end
    Face1 = [sommet1, sommet2];
    Face2 = [sommet2, sommet3];
    Face3 = [sommet1, sommet3];
    
    j = find_index(Face2,face);
    k = find_index(Face3,face);
    
    XF1 = XXF(i,:);
    XF2 = XXF(j,:);
    XF3 = XXF(k,:);
    
    X1=[Coorneu(sommet1,1),Coorneu(sommet1,2)];   %Coordonnées des sommets du triangle i
    X2=[Coorneu(sommet2,1),Coorneu(sommet2,2)];
    X3=[Coorneu(sommet3,1),Coorneu(sommet3,2)];
    
    X = [XF1 - Coorneu(sommet3,:)];
    n1 = [ - (X2(2) - X1(2)), X2(1) - X1(1)];
    n1 = n1/norm(n1);
    if (scalaire(n1,X)<0)
        n1 = - n1;
    end
    
    X = [XF2 - Coorneu(sommet1,:)];
    n2 = [ - (X3(2) - X2(2)), X3(1) - X2(1)];
    n2 = n2/norm(n2);
    if (scalaire(n2,X)<0)
        n2 = - n2;
    end
    
    X = [XF3 - Coorneu(sommet2,:)];
    n3 = [ -( X3(2) - X1(2)), X3(1) - X1(1)];
    n3 = n3/norm(n3);
    if (scalaire(n3,X) < 0)
        n3 = -n3;
    end
    
    K=abs(det([X2-X1;X3-X1]))/2;    %Mesure de K
    F1 = norm(X2 - X1);
    F2 = norm(X3 - X2);
    F3 = norm(X3 - X1);
    
    xk=[Xk(l,1),Xk(l,2)];               %barycentre
    
    
    Sgh1 = (Sgh(sommet1) + Sgh(sommet2))/2; %On prend Sh(xF) en moyennant sur les sommets adjacents de la face
    Sgh2 = (Sgh(sommet2) + Sgh(sommet3))/2;
    Sgh3 = (Sgh(sommet3) + Sgh(sommet1))/2;
    
    Ugh1 = Ugh(i);
    Ugh2 = Ugh(j);
    Ugh3 = Ugh(k);
    
   % Calcul du premier terme de l'estimateur d'erreur
    Y = 0; %car f = 0
    
    %Calcul du second terme de l'estimateur d'erreur
    
    Y = K * norm((Ugh1 - Sgh1)*F1/K*n1 + (Ugh2 - Sgh2)*F2/K*n2 + (Ugh3 - Sgh3)*F3/K*n3)^2;
    
    ETA_k2_2(l)=Y;   
end

eta_k2 = 0;
eta_k2_2 = 0;
for i = 1 : Nbtri
    eta_k2 = eta_k2 + ETA_k2(i,1);
    eta_k2_2 = eta_k2_2 + ETA_k2_2(i,1);
end
eta_k2 = sqrt(eta_k2);
eta_k2_2 = sqrt(eta_k2_2);


%%
%Plot des estimateurs d'erreurs

%Pour le cas 1 :
eta1 = [1.3829, 0.7408, 0.3976, 0.2202, 0.1283, 0.0791];
Vh1 = [20, 88, 368, 1504, 6080, 24448];

%Pour le cas 2
eta2 = [0.5551, 0.3256, 0.2394, 0.1699, 0.1189, 0.0829, 0.0577];
Vh1 = [14, 64, 272, 1120, 4544, 18304, 73472];


%%

%affiche2(UUh,Numtri,Coorneu,Numface,XXF,"Résolution dans le domaine [0,1]x[0,1] avec méthode non conforme");
%affiche(Sh,Numtri,Coorneu,"Résolution dans le domaine [0,1]x[0,1] avec methode non conforme");
%affiche(UU,Numtri,Coorneu);

%affiche2(Ugh,Numtri,Coorneu,Numface,XXF,"Résolution avec conditions non homogène ");
%affiche(Sgh,Numtri,Coorneu,"Résolution dans le domaine [-1,1]x[-1,1]\[0,1]x[0,1]");
%affiche(GG,Numtri,Coorneu, "Solution exacte");

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2022