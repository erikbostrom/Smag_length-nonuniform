function z=gll(Lz,Nz)
j = 0:Nz-1;
% GLL points in vertical
z = cos(pi*j/(Nz-1));
z = Lz*0.5*(1-z);
end