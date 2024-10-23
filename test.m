for l = 1:Nbtri
    i = Numtri(l,1);
    j = Numtri(l,2);
    k = Numtri(l,3);
    
    x = [Coorneu(i,1), Coorneu(j,1), Coorneu(k,1)];
    y = [Coorneu(i,2), Coorneu(j,2), Coorneu(k,2)];
    
    fill(x,y,sqrt(Lu(i)));
    colorbar;
    hold on;
end
                                                             
    