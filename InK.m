function bool= InK(i,Numtri,Coorneu,x,y)
bool=false;
M=[x,y];
S1=Numtri(i,1);                     %Sommets du triangle i
S2=Numtri(i,2);                        
S3=Numtri(i,3);
X1=[Coorneu(S1,1),Coorneu(S1,2)];   %Coordonnées des sommets du triangle i
X2=[Coorneu(S2,1),Coorneu(S2,2)];
X3=[Coorneu(S3,1),Coorneu(S3,2)];

P1=X1+(X2-X1)*scalaire(X2-X1,M-X1)/(norm(X2-X1)^2); %Projeté du point sur S1S2
P2=X1+(X3-X1)*scalaire(X3-X1,M-X1)/(norm(X1-X3)^2); %Projeté du point sur S1S3
P3=X2+(X3-X2)*scalaire(X3-X2,M-X2)/(norm(X2-X3)^2); %Projeté du point sur S2S3

%Normales$
A=X2-X1;         %S1S2
B=X3-X2;         %S2S3
C=X1-X3;         %S3S1
N1=[A(2),-A(1)];
if(scalaire(C,N1)<0)
    N1=-N1;
end
N2=[C(2),-C(1)];
if(scalaire(-A,N2)<0)
    N2=-N2;
end
N3=[B(2),-B(1)];
if(scalaire(A,N3)<0)
    N3=-N3;
end
if(scalaire(P1-M,N1)>=0 && scalaire(P2-M,N2)>=0 && scalaire(P3-M,N3)>=0)
    bool=true;
end
