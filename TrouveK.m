
function i = TrouveK(x,y,Nbtri,Numtri,Coorneu)
for k=1:Nbtri
    if (InK(k,Numtri,Coorneu,x,y)==true)
        i=k;
    end
end


        
    
