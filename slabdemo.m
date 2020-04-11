% INITIALIZE MATLAB
close all;
clc;
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FREE SPACE WAVELENGTH
lam0 = 1.0;
% SLAB PARAMETERS
n1 = 2.0;
n2 = 1.0;
a = 3*lam0;
% GRID
b = 5*lam0;
NRES = 10;
dx = lam0/NRES;
% NUMBER OF MODES TO CALCULATE
M = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD WAVEGUIDE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE GRID
Sx = a + 2*b;
Nx = ceil(Sx/dx);
Sx = Nx*dx;
%xa = [0.5:Nx-0.5]*dx;
%xa = xa - mean(xa);
% COMPUTE START AND STOP INDICES
nx = round(a/dx);
nx1 = round((Nx - nx)/2);
nx2 = nx1 + nx - 1;
% BUILD N
N = zeros(Nx,1);
N(1:nx1-1) = n2;
N(nx1:nx2) = n1;
N(nx2+1:Nx) = n2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM FD ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE k0
k0 = 2*pi/lam0;
% BUILD DX2
DX2 = sparse(Nx,Nx);
DX2 = spdiags(+1*ones(Nx,1),-1,DX2);
DX2 = spdiags(-2*ones(Nx,1), 0,DX2);
DX2 = spdiags(+1*ones(Nx,1),+1,DX2);
DX2 = DX2/ (k0*dx)^2;
% MAKE N DIAGONAL
N = diag(sparse(N(:)));
% SOLVE EIGEN-VALUE PROBLEM
A = DX2 + N^2;
[V,D] = eig(full(A));
NEFF = sqrt(diag(D));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VISUALIZE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SORT MODES
[~,ind] = sort(real(NEFF),'descend');
V = V(:,ind);
NEFF = NEFF(ind);
% OPEN FIGURE WINDOW
figure('Color','w');
hold on;
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
y = linspace(-b-a/2,b+a/2,Nx);
line(x,y,'Color','w','LineWidth',4);
h = line(x,y,'Color','b','LineWidth',2);
text(x0,y0,['Mode ' num2str(m)],'Rotation',-90,...
'HorizontalAlignment','center','VerticalAlignment','bottom');
text(x0,-y0,['n_{eff} = ' num2str(NEFF(m),'%4.2f')],'Rotation',-90,...
'HorizontalAlignment','center','VerticalAlignment','bottom');
end
% SET GRAPHICS VIEW
hold off;
h2 = get(h,'Parent');
set(h2,'XTick',[]);
axis equal tight;
ylabel('x axis');