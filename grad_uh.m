%%%%%%%% Implémentation de graduh %%%%%%%%%%%%%%%%%%%
function val = grad_uh(x,Nbtri,Numtri,Coorneu,face,UUh)
i=TrouveK(x(1),x(2),Nbtri,Numtri,Coorneu);


%Coordonnées des sommets

S1=Numtri(i,1);                     %Sommets du triangle i
S2=Numtri(i,2);                        
S3=Numtri(i,3);
X1=[Coorneu(S1,1),Coorneu(S1,2)];   %Coordonnées des sommets du triangle i
X2=[Coorneu(S2,1),Coorneu(S2,2)];
X3=[Coorneu(S3,1),Coorneu(S3,2)];

ja=find_index([S1,S2],face);        %indice de la face S1S2
jb=find_index([S2,S3],face);        %indice de la face S2S3
jc=find_index([S1,S3],face);        %indice de la face S3S3



A=X2-X1;         %S1S2
B=X3-X2;         %S2S3
C=X1-X3;         %S3S1

%Normales aux faces%
NA=[A(2),-A(1)];             %Normale à S1S2
if(scalaire(C,NA)<0)
    NA=-NA;
end
NA=NA/norm(NA);

NC=[C(2),-C(1)];             %Normale à S1S3
if(scalaire(-A,NC)<0)
    NC=-NC;
end
NC=NC/norm(NC);

NB=[B(2),-B(1)];             %Normale à S2S3
if(scalaire(A,NB)<0)
    NB=-NB;
end
NB=NB/norm(NB);
        
%Mesures%
K=det([A;B])/2; %mesure de K
FA=norm(A);     %mesure de la face S1S2
FB=norm(B);     %mesure de la face S2S3
FC=norm(C);     %mesure de la face S3S1

val=(UUh(ja)*FA*NA+UUh(jb)*FB*NB+UUh(jc)*FB*NB)/K;
end

