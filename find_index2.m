function j = find_index2(elem, Liste)
 for i = 1:size(Liste,1)
     if (elem == Liste(i,1) || elem == Liste(i,2) || elem == Liste(i,3))
         j = i;
     end
 end
 