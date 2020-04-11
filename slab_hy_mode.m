%initialize 
clc 
close all 
clear

%% developed by issah ibrahim 

%% sparse matrix 
M=8;
N =100; 
A = zeros(N,N); 
for i = 1:N
    A(i, i) = -2; 
end 

for i = 1:N-1
    A(i, i+1) = 1; 
    A(i+1, i) = 1; 
end 
da=.2;
lambda =1; 
ko = 2*pi/lambda; 
Dx2 = A/(ko*da)^2;

% dielectric 
Sx =1; 
eps = zeros(N, N); 
dx = Sx/N; 

a = .4; 
b = .4; 

nx = round(a/dx);
nx1 = 1+floor((N-nx)/2);
nx2 = nx1+nx -1; 
ny = round(b/dx); 
ny1 = 1+floor((N-ny)/2);
ny2 = ny1+ny -1; 

% refractive index 
n1 = 1; 
n2 = 2; 
eps(nx1:nx2, ny1:ny2) = 1; 

epsa = eps*n2 + (1-eps)*n1; 

epsa_act = diag(epsa); 
eps = epsa_act.*eye(N, N);


Amat = Dx2+eps^2; 

[V, D] = eig(Amat); 
NEFF = sqrt(diag(D)); 
[~,ind] = sort(real(NEFF),'descend');
V = V(:,ind);
NEFF = NEFF(ind);

a = 4; 
b = 4;
% DRAW SLAB WAVEGUIDE
x = [0 2*(M+1) 2*(M+1) 0 0];
y = [ -b-a/2 -b-a/2 b+a/2 b+a/2 -b-a/2 ];
fill(x,y,0.9*[1 1 1]);
y = [ -a/2 -a/2 a/2 a/2 -a/2 ];
fill(x,y,0.5*[1 1 1]);
% DRAW AND LABEL MODES
 
for m = 1 : M
x0 = 2*m;
y0 = (a + b)/2;
x = x0 + 3*V(:,m);
y = linspace(-b-a/2,b+a/2,N);
line(x,y,'Color','w','LineWidth',4);
h = line(x,y,'Color','b','LineWidth',2);
text(x0,y0,['Mode ' num2str(m)],'Rotation',-90,...
'HorizontalAlignment','center','VerticalAlignment','bottom');
text(x0,-y0,['n_{eff} = ' num2str(NEFF(m+1,1),'%4.2f')],'Rotation',-90,...
'HorizontalAlignment','center','VerticalAlignment','bottom');
end
% SET GRAPHICS VIEW
hold off;
h2 = get(h,'Parent');
set(h2,'XTick',[]);
axis equal tight;
ylabel('x axis');