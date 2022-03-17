% -------------------------------------------------------------------------
% Applying DMD to the solitons in the Nonlinear Schrodinger Equation (NLS)
%
% ADD DESCRIPTION HERE
%
% Author: Jason J. Bramburger
% -------------------------------------------------------------------------

% Clean workspace
clear all; close all; clc

%% Generate synthetic data from NLS using spectral methods

% Space
L = 40; n = 512;
Y = linspace(-L/2,L/2,n+1); x = Y(1:n);
k = (2*pi/L)*[0:n/2-1 -n/2:-1]';

% Time
slices = 20;
t = linspace(0,2*pi,slices+1); dt = t(2) - t(1);

% Initial condition
u = 2*sech(x);
ut = fft(u);

% Simulating the NLS
[t, utsol] = ode45(@(t,y) nls_rhs(y,k),t,ut);
for j = 1:length(t)
   usol(j,:) = ifft(utsol(j,:)); % back to x-space
end

%% Create X and Y data matrices

X = usol(1:end-1,:)';
Y = usol(2:end,:)';

%% Compute DMD matrix Y = A*X by A = Y * (pseduo-inverse of X) 

A = Y*pinv(X);
[eV, D] = eig(A); % compute eigenvalues + eigenvectors
mu = diag(D); % extract eigenvalues

% Continuous time eigenvalues for reference
omega = log(mu)/dt;

%% Plotting Discrete Eigenvalues (mu)

% make axis lines
line = -15:15;

% Unit circle
th = linspace(0,2*pi,1000);
xcos = cos(th);
ysin = sin(th);

figure(1)
plot(zeros(length(line),1),line,'k','Linewidth',2) % imaginary axis
hold on
plot(line,zeros(length(line),1),'k','Linewidth',2) % real axis
plot(xcos,ysin,'k--','LineWidth',2) % unit circle
plot(mu,'r.','Markersize',20)
xlabel('Re(\mu)')
ylabel('Im(\mu)')
set(gca,'FontSize',16,'Xlim',[-1.2 1.2],'Ylim',[-1.2 1.2])

%% Enforce DMD matrix is a unitary matrix
% --> Obtained through solution to Procrustes problem

% Solution to procrustes problem
M = X*Y';
[U_1, S_1, V_1] = svd(M, 'econ');
Au = U_1*V_1'; % unitary DMD matrix 

% Compute and plot eigenvalues
[eVu, Du] = eig(Au); % compute eigenvalues + eigenvectors
muu = diag(Du); % extract eigenvalues

figure(2)
plot(zeros(length(line),1),line,'k','Linewidth',2) % imaginary axis
hold on
plot(line,zeros(length(line),1),'k','Linewidth',2) % real axis
plot(xcos,ysin,'k--','LineWidth',2) % unit circle
plot(muu,'r.','Markersize',20)
xlabel('Re(\mu)')
ylabel('Im(\mu)')
set(gca,'FontSize',16,'Xlim',[-1.2 1.2],'Ylim',[-1.2 1.2])

%% Compare Forecasts far into the future

t = 0:dt:80;

% Regular DMD prediction
u_dmd = zeros(length(x),length(t));
u_dmd(:,1) = X(:,1);
for iter = 2:length(t)
   u_dmd(:,iter) = A*u_dmd(:,iter-1); 
end

% Unitary DMD matrix prediction
u_dmdu = zeros(length(x),length(t));
u_dmdu(:,1) = X(:,1);
for iter = 2:length(t)
   u_dmdu(:,iter) = Au*u_dmdu(:,iter-1); 
end

% Simulating the NLS far into the future
[t, utsollong] = ode45(@(t,y) nls_rhs(y,k),t,ut);
for j = 1:length(t)
   usollong(j,:) = ifft(utsollong(j,:)); % back to x-space
end

%% Plot prediction error

figure(3)
[X, T] = meshgrid(x,t);
surf(X,T,abs(u_dmdu' - usollong))
view(0,90)
shading interp
colorbar
set(gca,'FontSize',16,'Xlim',[-2 2])

%% NLS Right-Hand-Side

function rhs = nls_rhs(ut,k)
    u = ifft(ut);
    rhs = -(1i/2)*(k.^2).*ut + 1i*fft( (abs(u).^2).*u );
end
