function val = sigma_h(x,Nbtri,Numtri,Coorneu,face,UUh,fh,Xk)
val=-grad_uh(x,Nbtri,Numtri,Coorneu,face,UUh)+Fh(x,Nbtri,Numtri,Coorneu,fh,Xk);
end