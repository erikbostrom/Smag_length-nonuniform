% Compute the correction parameter gamma from
% Meyers & Sagaut 2006, 2007 for the SIMSON grid.
%
% requires the file gll.m
%
function comp_gamma

Nx = 48; Ny = 48; Nz = 33;
Lx = 2.5*pi; Ly = pi; Lz = 2;
[gamma,gammab,x,y,z]=f(Nx,Ny,Nz,Lx,Ly,Lz)
figure(1); plot(z,gammab,'k.-'); xlabel('z'); ylabel('\gamma'); hold on
ylim([0 1.2])

Nx = 48; Ny = 48; Nz = 49;
Lx = 2.5*pi; Ly = pi; Lz = 2;
[gamma,gammab,x,y,z]=f(Nx,Ny,Nz,Lx,Ly,Lz)
figure(1); plot(z,gammab,'k.-.'); xlabel('z'); ylabel('\gamma');
ylim([0 1.2])

Nx = 48; Ny = 48; Nz = 65;
Lx = 2.5*pi; Ly = pi; Lz = 2;
[gamma,gammab,x,y,z]=f(Nx,Ny,Nz,Lx,Ly,Lz)
figure(1); plot(z,gammab,'k.--'); xlabel('z'); ylabel('\gamma');
ylim([0 1.2])

Nx = 48; Ny = 48; Nz = 97;
Lx = 2.5*pi; Ly = pi; Lz = 2;
[gamma,gammab,x,y,z]=f(Nx,Ny,Nz,Lx,Ly,Lz)
figure(1); plot(z,gammab,'b.-'); xlabel('z'); ylabel('\gamma');
ylim([0 1.2])

hold off

end

function [gamma,gammab,x,y,z]=f(Nx,Ny,Nz,Lx,Ly,Lz)
x  = linspace(0,Lx,Nx);
y  = linspace(0,Ly,Ny);
z  = gll(Lz,Nz);
dx = Lx/Nx;
dy = Ly/Ny;
dz = [z(2) sqrt((z(3:end)-z(2:end-1)).*(z(2:end-1)-z(1:end-2))) z(2)];
delta = (dx*dy.*dz).^(1/3);

xi = linspace(-1,1,48);
a  = 1;
b  = 1;
zz = a*tanh(b*xi);

%z2 = gll(Lz,Nz+2);
%dz = sqrt((z2(3:end)-z2(2:end-1)).*(z2(2:end-1)-z2(1:end-2)))

%figure; plot(z,dz,'.-'); xlabel('z'); ylabel('h_z');
disp('The sum of hz should be close to two: '); disp(sum(dz));
disp('In the Meyers & Sagaut 2007 paper they have hy0=1/500=0.02 and hyc=0.049.'); disp(dz(1)); disp(dz(round(end/2)));

kx = pi/dx;
ky = pi/dy;
kz = pi./dz

% Integrate ellipse
for i=1:length(kz)
f = @(theta,phi) (...
    (cos(theta).^2.*sin(phi).^2)./(kx^2) + ...
    (sin(theta).^2.*sin(phi).^2)./(ky^2) + ...
    (cos(phi).^2)./(kz(i)^2)).^(-2/3)...
    .*sin(theta);
a(i) = integral2(f,0,pi,0,2*pi)
end
gamma = (1/(3*pi)*a./(pi./delta).^(4/3)).^(3/4);
disp(a);
%figure; plot(z,a)
%figure; plot(z,gamma)

% Integrate box
g = @(x,y,z) (x.^2 + y.^2 + z.^2).^(-5/6);
for i = 1:length(kz)
b1 = integral3(g,0,kx,0,ky,0,kz(i));
b2(i) = 8*b1;
end
gammab = (1/(3*pi)*b2./((pi./delta).^(4/3))).^(3/4);
%integral3(g,0,kx,0,ky,0,kz(1))
%integral3(g,0,kx,0,ky,-kz(1),0)
%integral3(g,0,kx,-ky,0,0,kz(1))
%integral3(g,0,kx,-ky,0,-kz(1),0)
%integral3(g,-kx,0,0,ky,0,kz(1))
%integral3(g,-kx,0,0,ky,-kz(1),0)
%integral3(g,-kx,0,-ky,0,0,kz(1))
%integral3(g,-kx,0,-ky,0,-kz(1),0)
% All the same...
end
%b = 3*pi*(pi./delta).^(4/3);