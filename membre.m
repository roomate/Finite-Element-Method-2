function bool = membre(elem, liste) %Regarde si un élément fait bien partie de la liste
bool = false;
for i = 1:size(liste,1)
    if ((elem(1) == liste(i,1) && elem(2) == liste(i,2)) || (elem(2) == liste(i,1) && elem(1) == liste(i,2)))
        bool = true;
    end
end
