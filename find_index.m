function index = find_index(elem,face) %%Trouver le vecteur elem 1x2 dans la liste Face Nbface_int*2 
index = 1;
if (membre(elem,face) == 0)
    disp("oulah");
end
for j=1:size(face,1)
    if ((face(j,1) == elem(1) && face(j,2) == elem(2)) || (face(j,2) == elem(1) && face(j,1) == elem(2)))
        index=j;
    end
end    
    