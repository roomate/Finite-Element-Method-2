function val = Fh(x,Nbtri,Numtri,Coorneu,fh,Xk)
i=TrouveK(x(1),x(2),Nbtri,Numtri,Coorneu);         % on se place dans la bonne maille
val=(fh(i,1)/2)*([x(1)-Xk(i,1),x(2)-Xk(i,2)]);
end

