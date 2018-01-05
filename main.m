clear all;
hold all;

Nx = 48;
Ny = 48;
Nz = 97;

Lx = 2*pi;
Ly = pi;
Lz = 2;

x  = linspace(0,Lx,Nx);
y  = linspace(0,Ly,Ny);
z  = gll(Lz,Nz);
dx = Lx/Nx;
dy = Ly/Ny;
dz = [z(2) sqrt((z(3:end)-z(2:end-1)).*(z(2:end-1)-z(1:end-2))) z(2)];

delta = (dx*dy.*dz).^(1/3);

kx = pi/dx;
ky = pi/dy;
kz = pi./dz;

ky = kx;
kz = ones(length(kz),1)*kx;
delta = dx;

% Integrate box
g = @(x,y,z) (x.^2 + y.^2 + z.^2).^(-5/6);
for i = 1:length(kz)
    b1 = integral3(g,0,kx,0,ky,0,kz(i));
    b2(i) = 8*b1;
end
gammab =(1/(3*pi)*b2./((pi./delta).^(4/3))).^(3/4);

%% Plotting
figure(1); plot(z,gammab); ylim([0 1.3]);