% -------------------------------------------------------------------------
% Proper Orthogonal Decomposition of a Spiral Wave 
%
% This code applies proper orthogonal decomposition (POD) to spacet-time 
% data obtained from numerically integrating a reaction-diffusion PDE. In
% particular we focus on a spiral wave, which represents a spatially 
% coherent, time-periodic structure that simply rotates about a fixed point
% in space. The POD decomposition shows that most of the `energy' is 
% contained in the first two modes. These modes can be used to faithfully
% reconstruct the space-time dynamics of the spiral at a significant
% reduction in the amount of data stored.
%
% This script accompanies Section 1.# of Data-Driven Methods for
% Dynamic Systems. 
%
% Author: Jason J. Bramburger
% -------------------------------------------------------------------------

% Clean workspace
clear all; close all; clc

%% General spiral wave data using spectral methods

% System parameters
d1 = 0.1; 
d2 = 0.1; 
beta = 1.0;

% Spectral method variables
tspan = 0:0.05:40;
L = 20; 
n = 256; 
N = n*n;
x2 = linspace(-L/2,L/2,n+1); 
x = x2(1:n); 
y = x;
kx = (2*pi/L)*[0:(n/2-1) -n/2:-1]; 
ky = kx;

% Initial Conditions
m = 1; % number of spirals
[X,Y] = meshgrid(x,y);
[KX,KY] = meshgrid(kx,ky);
K2 = KX.^2 + KY.^2; 
K22=reshape(K2,N,1);

u = zeros(length(x),length(y),length(tspan));
v = zeros(length(x),length(y),length(tspan));

u0 = tanh(sqrt(X.^2+Y.^2)).*cos(m*angle(X+1i*Y)-(sqrt(X.^2+Y.^2)));
v0 = tanh(sqrt(X.^2+Y.^2)).*sin(m*angle(X+1i*Y)-(sqrt(X.^2+Y.^2)));

% Numerical timestepping in Fourier domain
uv0 = [reshape(fft2(u0),1,N) reshape(fft2(v0),1,N)].';
[t, uvsol] = ode45(@(t,uvt) rdpde(t,uvt,K22,d1,d2,beta,n,N),tspan,uv0);

%% Reshape and view spiral wave

for j = 1:length(tspan)
    ut = reshape((uvsol(j,1:N).'),n,n);
    vt = reshape((uvsol(j,(N+1):(2*N)).'),n,n);
    u(:,:,j) = real(ifft2(ut));
    v(:,:,j) = real(ifft2(vt));

    figure(1)
    pcolor(x,y,u(:,:,j)); shading interp; colormap(hot); colorbar; drawnow; 
end

%% Apply SVD to space-time data

% Aggregate data
u = reshape(u,n^2,length(tspan));
v = reshape(v,n^2,length(tspan));
F = [u; v]; % data matrix

% Apply SVD
[U, S, V] = svd(F,'econ');

%% Plot singular values

figure(2)
plot(diag(S.^2)/sum(diag(S.^2)),'ko','MarkerSize',10)
xlabel('$k$','interpreter','latex','FontWeight','Bold')
ylabel('$\sigma_k^2/\sum\sigma_j^2$','interpreter','latex','FontWeight','Bold')
set(gca,'Fontsize',16,'Xlim',[0.9 50])

%% Plot POD modes

% First three POD modes
pod1 = reshape(U(1:n^2,1),n,n);
pod2 = reshape(U(1:n^2,2),n,n);
pod3 = reshape(U(1:n^2,3),n,n);

% Plot mode 1
figure(3)
surf(X,Y,pod1)
view(0,90)
shading interp
xlabel('$x$','interpreter','latex','FontWeight','Bold')
ylabel('$y$','interpreter','latex','FontWeight','Bold')
set(gca,'Fontsize',16)

% Plot mode 2
figure(4)
surf(X,Y,pod2)
view(0,90)
shading interp
xlabel('$x$','interpreter','latex','FontWeight','Bold')
ylabel('$y$','interpreter','latex','FontWeight','Bold')
set(gca,'Fontsize',16)

% Plot mode 3
figure(5)
surf(X,Y,pod3)
view(0,90)
shading interp
xlabel('$x$','interpreter','latex','FontWeight','Bold')
ylabel('$y$','interpreter','latex','FontWeight','Bold')
set(gca,'Fontsize',16)

%% Plot time evolution of POD modes

% First 3 temporal modes
a1 = V(:,1);
a2 = V(:,2);
a3 = V(:,3);

% Plot modes
figure(6)
plot(tspan,a1,'Color',[1 69/255 79/255],'LineWidth',2)
hold on
plot(tspan,a2,'Color',[36/255 122/255 254/255],'LineWidth',2)
plot(tspan,a3,'Color',[0 200/255 140/255],'LineWidth',2)
xlabel('$t$','interpreter','latex','FontWeight','Bold')
legend('$V_1(t)$','$V_2(t)$','$V_3(x)$','Location','Best','interpreter','latex')
set(gca,'Fontsize',16,'Xlim',[0 40])

%% Low-rank reconstructions

rank1 = S(1,1)*U(:,1)*V(end,1);
rank2 = rank1 + S(2,2)*U(:,2)*V(end,2);
rank3 = rank2 + S(3,3)*U(:,3)*V(end,3);
rank1 = reshape(rank1(1:n^2),n,n);
rank2 = reshape(rank2(1:n^2,1),n,n);
rank3 = reshape(rank3(1:n^2,1),n,n);

% Plot final panel of the true solution
figure(7)
surf(X,Y,reshape(u(:,end),n,n))
view(0,90)
shading interp
xlabel('$x$','interpreter','latex','FontWeight','Bold')
ylabel('$y$','interpreter','latex','FontWeight','Bold')
set(gca,'Fontsize',16)

% Rank 1 Reconstruction
figure(8)
surf(X,Y,rank1)
view(0,90)
shading interp
xlabel('$x$','interpreter','latex','FontWeight','Bold')
ylabel('$y$','interpreter','latex','FontWeight','Bold')
set(gca,'Fontsize',16)

% Rank 2 Reconstruction
figure(9)
surf(X,Y,rank2)
view(0,90)
shading interp
xlabel('$x$','interpreter','latex','FontWeight','Bold')
ylabel('$y$','interpreter','latex','FontWeight','Bold')
set(gca,'Fontsize',16)

% Rank 3 Reconstruction
figure(10)
surf(X,Y,rank3)
view(0,90)
shading interp
xlabel('$x$','interpreter','latex','FontWeight','Bold')
ylabel('$y$','interpreter','latex','FontWeight','Bold')
set(gca,'Fontsize',16)

%% Reaction-diffusion RHS

function rhs = rdpde(t,uvt,K22,d1,d2,beta,n,N)

    % Calculate u and v terms
    ut = reshape((uvt(1:N)),n,n);
    vt = reshape((uvt((N+1):(2*N))),n,n);
    u = real(ifft2(ut)); 
    v = real(ifft2(vt));

    % Reaction Terms
    u3 = u.^3; 
    v3 = v.^3; 
    u2v = (u.^2).*v; 
    uv2 = u.*(v.^2);
    utrhs = reshape((fft2(u-u3-uv2+beta*u2v+beta*v3)),N,1);
    vtrhs = reshape((fft2(v-u2v-v3-beta*u3-beta*uv2)),N,1);

    rhs = [-d1*K22.*uvt(1:N)+utrhs
         -d2*K22.*uvt(N+1:end)+vtrhs];
 
end