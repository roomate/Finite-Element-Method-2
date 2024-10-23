%%%%%%%%%%%%%%%%%%%%% EXO 4 %%%%%%%%%%%%%%%%%%%%%%%%%%

ETA_k2=zeros(Nbtri,1);
Sh = zeros();
%% 
for i=1:Nbtri
    S1=Numtri(i,1);       %sommets du triangle i
    S2=Numtri(i,2); 
    S3=Numtri(i,3);
    
    X1=[Coorneu(S1,1),Coorneu(S1,2)];   %Coordonnées des sommets du triangle i
    X2=[Coorneu(S2,1),Coorneu(S2,2)];
    X3=[Coorneu(S3,1),Coorneu(S3,2)];
    
    K=det([X2-X1;X3-X2])/2;             %Mesure de K
    
    M1=(X1+X2)/2;                       %Milieux des côtés                   
    M2=(X1+X3)/2;
    M3=(X2+X3)/2;
    xk=[Xk(i,1),Xk(i,2)];               %barycentre
    
    Y = ((K/3)*fh(i,1)^2/4)*( scalaire(M1-xk,M1-xk) + scalaire(M2-xk,M2-xk) + scalaire(M3-xk,M3-xk));
    
    ETA_k2(i,1)=Y;
   
    
end

eta_k2=0;
for i=1:Nbtri
    eta_k2=eta_k2+ETA_k2(i,1);
end