% -------------------------------------------------------------------------
% Identifying Optimal Controllers from Data
%
% The goal of this script is to use data-driven approximations of the Lie
% derivative to identify optimal controllers from data. This script uses
% the inverted pendulum on a cart model
% 
%         theta'' = sin(theta) - eps*theta' - cos(theta)u
%
% where eps > 0 is the friction coefficient and u is the control input. 
% Here the control is given by the acceleration of the cart. 
%
% This script identifies a state-dependent feedback controller 
% u = u_*(theta,theta') that stabilizes the pendulum in the upright
% position, i.e. at the equilibrium (theta,theta') = (0,0). The control law
% u_* is identified through polynomial optimization with a control Lyapunov
% function. 
%
% Application of the method transforms the state variables theta and theta'
% to the 3D observables: x_1 = cos(theta), x_2 = sin(theta), and x_3 =
% theta'. We then discover u_* as a polynomial function of (x_1,x_2,x_3).
% 
% To run this notebook one requires paths to the freely available software
% packages YALMIP and MOSEK. YALMIP can be downloaded at:
%           https://yalmip.github.io/download/
% and MOSEK can be downloaded at:
%           https://www.mosek.com/downloads/
%
% This script accompanies Section 4.4 of Data-Driven Methods for
% Dynamic Systems.
%
% Written by J. Bramburger 
%
% -------------------------------------------------------------------------

% Clean workspace
clear; close all; clc
yalmip clear
format long

%% Method Parameters 
% maxDp = max degree of Dp dictionary of obserables
% maxDq = max degree of Dq dictionary of obserables
% degU = degree of the controller u as a function of x variables
% al = decay parameter for Lyapunov function
% eta = state space parameter
maxDp = 3;
maxDq = 4;
degU = 4;
al = 100;
eta = sqrt(0.95);

%% Generated Controlled Data
dt = 0.01;
t = 0:dt:20;

% Function rhs - we use random sinusoidal control inputs to learn the
%        effect of the control on the data
syst = @(t,x,a,b) [x(2); -0.1*x(2) + sin(x(1)) - cos(x(1))*a*sin(t + b)];

% Initialize for multiple initial conditions
numICs = 20;
xdat = [];
ydat = [];
uxdat = [];
uydat = [];
for ind = 1:numICs

    % Randomized initial conditions
    x0 = [2*pi*rand - pi; 4*rand - 2];
    
    % ODE45 to simulate ODE
    a = 2*rand-1;
    b = 2*pi*rand - pi;
    rhs = @(t,x) syst(t,x,a,b);
    [t, sol] = ode45(rhs,t,x0);
    
    sol = [cos(sol(:,1)), sin(sol(:,1)), sol(:,2)];
    
    % Aggregate data
    xdat = [xdat, sol(1:end-1,:)'];
    ydat = [ydat, sol(2:end,:)'];
    
    % Control 
    uxdat = [uxdat; a*sin(t(1:end-1)+b)];
    uydat = [uydat; a*sin(t(2:end)+b)];
end

% Transpose
xdat = xdat';
ydat = ydat';

%% Create P and Qu matrices

% Qu matrix
pow = monpowers(3,maxDq);
m = size(pow,1); % number of nontrivial monomials in Q
Q = [];
for i = 1:m
   
    zx = xdat.^pow(i,:); 
    
    if pow(i,2) >= 2 % no polynomials with deg(x_2) > 2
        Q(i,:) = zeros(1,length(zx));    
    else
        Q(i,:) = prod(zx,2);
    end
end

% Append in control to make Qu matrix
Q = [Q; Q.*uxdat'];

% P matrix
pow = monpowers(3,maxDp);
ell = size(pow,1); % number of nontrivial monomials in P
P = [];
for i = 1:ell
    
    zy = ydat.^pow(i,:);
   
    if pow(i,2) >= 2 % no polynomials with deg(x_2) > 2  
        P(i,:) = zeros(1,length(zy));
    else
        P(i,:) = prod(zy,2)';
    end
end

%% Create Koopman and Lie derivative matrix approximations

% Koopman approximation
K = P*pinv(Q);

% Symbolic variables
x = sdpvar(3,1); % 2D state variable
z = monolist(x,maxDp,0); % monomials that make up span(Dp)
w = monolist(x,maxDq,0); % monomials that make up span(Dq)

% Controller
[u,cu] = polynomial(x,degU);

% Lyapunov function vector coefficients in span(Dp)
c = zeros(length(z),1);
c(1) = 1 + al;
c(2) = -1;
c(10) = 0.5;
c(11) = -al;

% Lyapunov function
% Lyap = 0.5*( 0.5*x(3)^2 + x(1) - 1 )^2 + al*(1 - x(1)^3);
% c = coefficients(Lyap,x);

% Lie derivative approximation
L = (K - eye( size(K) ) ) /dt;
thresh = 0.05; % use to stamp out noise
L(abs(L) <= thresh) = 0;
Lie = c.'*(L(:,1:m)*w + L(:,m+1:end)*w*u);

%% Polynomial optimization to find controller 

% S procedure
[s1,cs1] = polynomial(x,degU,1);
[s2,cs2] = polynomial(x,degU,1);

% Lyapunov function used:
%      V = [0.5*th'^2 + 1 - cos(th)] + al*(1 - cos(th)^3) 
obj = sos( -Lie + (1 - x(1)^2 - x(2)^2)*s1 - (eta^2 - x(2)^2)*s2); 
obj = [obj; sos(s2)];

% SOS solver
ops = sdpsettings('solver','mosek','verbose',1); % verbose = 1/0 --> YALMIP display on/off
sol = solvesos(obj,sum(abs(cu)),ops,[cu; cs1; cs2])
flag = sol.problem;

%% Clean small values and print controller

% Remove coefficients smaller than 1e-4
cu = clean(value(cu),1e-4);

% Print controller
xmons = monolist(x,degU);
fprintf('\nIdentified control input: \n')
fprintf('u(x) = ')
sdisplay(value(cu).'*xmons)
fprintf('\n')

%% Check that we have indeed discovered a control law

% Randomized initial conditions
x0 = [2*pi*rand - pi; 2*rand - 1];

% Initial condition close to hanging down position
x0 = [pi - 0.001; 0];

% Integrate controlled system
upow = monpowers(3,degU);
cuv = value(cu);
dt = 0.01;
t = 0:dt:8;
rhs = @(t,x) systc(t,x,cuv,upow);
[t,xc] = ode45(rhs,t,x0);

% Plot controlled solution
clf
figure(1)
plot(t,xc(:,1),'Color',[36/255 122/255 254/255],'LineWidth',3)
hold on
plot(t,xc(:,2),'--','Color',[1 69/255 79/255],'LineWidth',3)
xlabel('$t$','interpreter','latex')
title('Controlled Solutions')
legend('$\theta(t)$','$\dot\theta(t)$','Location','Best','interpreter','latex','FontSize',20)
set(gca,'FontSize',16,'Xlim',[0,t(end)])
axis tight

% Determine controller as time goes on
uc = [];
for j = 1:length(t)
   uc(j,1) = cuv.'*prod([cos(xc(j,1)) sin(xc(j,1)) xc(j,2)].^upow,2);  
end

figure(2)
V0 = 10; % (0.5*xc(1,2).^2 + 1 - cos(xc(1,1)) + al*(1 - cos(xc(1,1)).^3));
uc0 = 1; % uc(1);
plot(t,(0.5*xc(:,2).^2 + 1 - cos(xc(:,1)) + al*(1 - cos(xc(:,1)).^3))/V0,'Color',[36/255 122/255 254/255],'LineWidth',3)
hold on
plot(t,uc/uc0,'--','Color',[1 69/255 79/255],'LineWidth',3)
xlabel('$t$','interpreter','latex')
title('Controller and Lyapunov Function')
legend('$V(\theta(t),\dot\theta(t))$','$u(t)$','Location','Best','interpreter','latex','FontSize',20)
set(gca,'FontSize',16,'Xlim',[0,t(end)])
axis tight



%% Controlled ODE system
function dxc = systc(t,x,cu,upow)

    dxc(1) = x(2);
    dxc(2) = -0.1*x(2) + sin(x(1)) - cos(x(1))*cu.'*prod([cos(x(1)) sin(x(1)) x(2)].^upow,2);  
    
    dxc = dxc';
end


