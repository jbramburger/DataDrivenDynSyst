% -------------------------------------------------------------------------
% Sliding Window DMD
%
% TBD.
%
% This script accompanies Section 2.2 of Data-Driven Methods for
% Dynamic Systems. 
%
% Author: Jason J. Bramburger
% -------------------------------------------------------------------------

% Clean workspace
clear all; close all; clc

%% Generate data 

% Model parameters
eps = 0.01;
delta = 10;

% ODE step parameters
dt = eps/100; % need lots of resolution for fast timescale
tspan = 0:dt:2^6;

% Initial conditions
M = 4; % dimension of ODE
x0 = [0.0; 0.5; 0.0; 0.5];

% Integrate ODE
[t, x] = ode45(@(t,x) rhs(x,delta,eps),tspan,x0);

%% Random mixing of the singals

% Generate random orthogonal matrix
Q = rand(M);
Q = orth( Q.' ).'; % orthogonal rows
Q = Q*diag(2*randi([0, 1],[1 M])-1); % randomly make some components negative as well

% To reproduce results from the textbook load orth_mat.mat
load orth_mat.mat

% Mix signal using the orthogonal matrix
xdat = x*Q;

% Plot multiscale signal
figure(1)
plot(t,xdat(:,1),'k','LineWidth',1)
hold on
plot(t,xdat(:,2),'k','LineWidth',1)
plot(t,xdat(:,3),'k','LineWidth',1)
plot(t,xdat(:,4),'k','LineWidth',1)
xlabel('$t$','interpreter','latex','FontSize',16)
ylabel('$x(t)$','interpreter','latex','FontSize',16)
set(gca,'FontSize',14,'Xlim',[0,tspan(end)])

%% Apply standard DMD to multiscale signal

% Create X and Y matrices
X = xdat(1:end-1,:)';
Y = xdat(2:end,:)';

% DMD matrix
A = Y*pinv(X);
[eV, D] = eig(A); % compute eigenvalues + eigenvectors
mu = diag(D); % extract eigenvalues

% Continuous time eigenvalues for reference
omega = log(mu)/dt;

%% Forecast with DMD

N = length(t);
xforecast = zeros(M,N);
xforecast(:,1) = X(:,1);

% DMD forecast
for m = 1:N-1
    xforecast(:,m+1) = A*xforecast(:,m);
end

% Compare DMD forecast with data
% --> Plotting x_3(t) because the difference is easiest to see
figure(2)
plot(t,xdat(:,3),'k','LineWidth',1)
hold on
plot(t,xforecast(3,:),'Color',[1 69/255 79/255],'LineWidth',1)
xlabel('$t$','interpreter','latex','FontSize',16)
ylabel('$x_3(t)$','interpreter','latex','FontSize',16)
set(gca,'FontSize',14,'Xlim',[0,tspan(end)])

%% ------------------------------------------------------- 
% Begin Windowed DMD
% -------------------------------------------------------

%% 



%% ODE right-hand-side

function dx = rhs(x,delta,eps)
    x1 = x(1); 
    x2 = x(2); 
    y1 = x(3); 
    y2 = x(4);
    
    dx=[x2; -y1^2 * x1^3; y2; -eps^(-1)*y1 - delta*y1^3];
end