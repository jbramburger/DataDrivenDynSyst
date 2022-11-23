% -------------------------------------------------------------------------
% Kernel Dynamic Mode Decomposition (EDMD)
%
% The goal of this script is to apply the kernel method of dynamic mode
% decomposition to solutions to the Burgers equation, given by the
% nonlinear partial differential equation:
%             u_t = nu*u_xx - u*u_x.
% The method is meant to circumvent exploding dictionary sizes coming from
% high-dimensional data by replacing inner products of dictionary
% evaluations with an evaluation of a kernel function. Here we use the
% polynomial kernel h(x,y) = (1 + y^T*x)^d, where we take the degree d to
% be a parameter. 
%
% The last block in the code applies EDMD to the Burgers
% equation using the single observable given by the Cole-Hopf
% transformation. The significance is that this transformation v(u) turns the
% Burgers equation into the linear heat equation
%             v_t = nu*v_xx
% and so the Koopman expansion can be full realized using separation of
% variables. 
%
% This script accompanies Section 2.5 of Data-Driven Methods for
% Dynamic Systems. 
%
% Author: Jason J. Bramburger
% -------------------------------------------------------------------------

% Clean workspace
clear all; close all; clc

%% Simulate Burgers equation

% Viscosity parameter
nu = 0.1;

% Spatial and temporal domain parameters
dt = 1e-1;
t = 0:dt:10; 
L = 2*pi; % spatial domain is [-L/2,L/2] - L is the length
n = 2^8; % spatial points in each dimension
x = linspace(-L/2,L/2,n);
h = x(2) - x(1);

% d_x with Dirichlet BCs
D = sparse(1:n-1,[2:n-1 n],ones(n-1,1)/2,n,n); 
D = (D - D')/h;

% d_xx with Dirichlet BCs
e = ones(n,1);
D2 = sparse(1:n-1,[2:n-1 n],ones(n-1,1),n,n) - sparse(1:n,[1:n],e,n,n);
D2 = (D2 + D2');
D2 = D2/h^2;

% Initialize for multiple initial conditions
numICs = 5;
xdat = [];
ydat = [];
for ind = 1:numICs

    % Randomized initial conditions
    numSines = randi([1,10]);
    u0 = zeros(1,n);
    for jnd = 1:numSines
       u0 = u0 + (2*rand - 1)*sin(jnd*x); 
    end
    
    % Simulate PDE using Fourier methods
    [t,usol] = ode45(@(t,u) burgers(t,u,nu,D,D2),t,u0);
    
    % Aggregate data
    xdat = [xdat, usol(1:end-1,:)'];
    ydat = [ydat, usol(2:end,:)'];
end

% Append zero solution
% xdat = [xdat, zeros(n,1)];
% ydat = [ydat, zeros(n,1)];

%% Plot solution for visualization

% Spatial domain
[X, T] = meshgrid(x,t);

figure(1)
hold on
surf(T',X',usol')
plot3(0*x,x,usol(1,:),'Color',[1 69/255 79/255],'LineWidth',5)
plot3(5*ones(1,n),x,usol(51,:),'Color',[36/255 122/255 254/255],'LineWidth',5)
hold off
grid on
view(-45,25)
shading interp
colormap gray
xlabel('t')
ylabel('x')
zlabel('u(x,t)')
set(gca,'FontSize',16,'Xlim',[0 5],'Ylim',[-L/2 L/2])

%% Build data matrices with kernel

% Number of time elements
N = length(xdat(1,:)); 

% Polynomial degree
d = 8;

% Kernel matrices
PsiXPsiX = (ones(N,N) + xdat.'*xdat).^d;
PsiYPsiX = (ones(N,N) + ydat.'*xdat).^d;

%% Diagonalize PsiXPsiX using SVD to get V and Sigma

r = N; % Decrease rank to speed up computations
[U, Sigma2, V] = svds(PsiXPsiX,r);
Sigma = real(sqrt(Sigma2));

%% Create hat(A) matrix and analyze spectrum

Ahat = pinv(Sigma)*U.'*PsiYPsiX*V*pinv(Sigma);
[what, mu, vhat] = eig(Ahat);

% make axis lines
line = -15:15;

% Unit circle
th = linspace(0,2*pi,1000);
xcos = cos(th);
ysin = sin(th);

% Plot kernel DMD eigenvalues
figure(2)
clf
plot(zeros(length(line),1),line,'k','Linewidth',2) % imaginary axis
hold on
plot(line,zeros(length(line),1),'k','Linewidth',2) % real axis
plot(xcos,ysin,'k--','LineWidth',2) % unit circle
plot(diag(mu),'.','Color',[1 69/255 79/255],'Markersize',30)
xlabel('Re(\mu)')
ylabel('Im(\mu)')
set(gca,'FontSize',16,'Xlim',[-1.2 1.2],'Ylim',[-1.2 1.2])

%% DMD on Cole-Hopf transformation

%  Create Cole-Hopf transformed variable v(u)
xi = zeros(length(t),n);
vCH = zeros(length(t),n);
for ind = 2:n-1
    xi(:,ind) = (-1/(2*nu))*trapz(x(1:ind),usol(:,1:ind),2)';
end
expxi = exp(xi);
expxiInt = trapz(x,expxi,2);
vCH = expxi./expxiInt;

% Create DMD matrices
X = vCH(1:end-1,:)';
Y = vCH(2:end,:)';

% DMD matrix
A = Y*pinv(X);
[eV, D] = eig(A); % compute eigenvalues + eigenvectors
mu = diag(D); % extract eigenvalues

% Continuous time eigenvalues for reference
omega = log(mu)*L^2/dt/nu/pi^2;

% Check which are close to integers
sort(sqrt(-omega))

%% Burgers right-hand-side
function rhs = burgers(t,u,nu,D,D2)
    
    % Create nonlinear term
    uux = u.*D*u;
    
    % Burgers equation in Fourier space
    rhs = nu*D2*u - uux;

end








