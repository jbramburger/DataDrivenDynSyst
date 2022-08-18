% -------------------------------------------------------------------------
% Discovering Conserved Quantities
% -------------------------------------------------------------------------
%
% The goal of this script is to discover conserved quantities from data
% gather from a dynamical system. We illustrate with a simple 3D dynamical
% system that has two conserved quantities: one whose level sets are
% spheres and the other whose are ellipsoids. The method using the SVD to 
% identify orthogonal basis vectors that span the nullspace of Gamma, as 
% defined in the text. 
%
% This script accompanies Section 3.5 of Data-Driven Methods for
% Dynamic Systems. 
%
% Author: Jason J. Bramburger
% -------------------------------------------------------------------------

% Clean workspace
clear all; close all; clc

%% Generate training data

% Integration parameters
dt = 0.01;
t = 0:dt:5;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,3));

% Simulate conservative system
x = [];
numInit = 20; % number of initial conditions
for n = 1:numInit
    x0 = 10*rand(1,3)-5; 
    [~,out] = ode45(@(t,x) consODE(x),t,x0,options);
    x = [x; out];
end

%% Obtain derivative estimates from ODE right-hand-side

for j = 1:length(x(:,1))
    dxdt(j,:) = consODE(x(j,:)); % ODE rhs
end

%% Create Gamma(X,dotX) matrix

Gamma(:,1:3) = dxdt; % x1, x2, & x3 observables
Gamma(:,4) = 2*x(:,1).*dxdt(:,1); % x1^2 observable
Gamma(:,5) = x(:,2).*dxdt(:,1) + x(:,1).*dxdt(:,2); % x1*x2 observable
Gamma(:,6) = x(:,3).*dxdt(:,1) + x(:,1).*dxdt(:,3); % x1*x3 observable
Gamma(:,7) = 2*x(:,2).*dxdt(:,2); % x2^2 observable
Gamma(:,8) = x(:,3).*dxdt(:,2) + x(:,2).*dxdt(:,3); % x2*x3 observable
Gamma(:,9) = 2*x(:,3).*dxdt(:,3); % x3^2 observable
Gamma(:,10) = 3*x(:,1).*x(:,1).*dxdt(:,1); % x1^3 observable
Gamma(:,11) = 2*x(:,1).*x(:,2).*dxdt(:,1) + x(:,1).*x(:,1).*dxdt(:,2); % x1*x1*x2 observable
Gamma(:,12) = x(:,2).*x(:,2).*dxdt(:,1) + 2*x(:,1).*x(:,2).*dxdt(:,2); % x1*x2*x2 observable
Gamma(:,13) = x(:,2).*x(:,3).*dxdt(:,1) + x(:,1).*x(:,3).*dxdt(:,2) + x(:,1).*x(:,2).*dxdt(:,3); % x1*x2*x3 observable
Gamma(:,14) = x(:,3).*x(:,3).*dxdt(:,1) + 2*x(:,1).*x(:,3).*dxdt(:,3); % x1*x3*x3 observable
Gamma(:,15) = 3*x(:,2).*x(:,2).*dxdt(:,2); % x2^3 observable
Gamma(:,16) = 2*x(:,2).*x(:,3).*dxdt(:,2) + x(:,2).*x(:,2).*dxdt(:,3); % x2*x2*x3 observable
Gamma(:,17) = x(:,3).*x(:,3).*dxdt(:,2) + 2*x(:,2).*x(:,3).*dxdt(:,3); % x2*x3*x3 observable
Gamma(:,18) = 3*x(:,3).*x(:,3).*dxdt(:,3); % x3^3 observable

%% Compute SVD of Gamma and analyze zero singular values

[U, S, V] = svd(Gamma);

% Plot singular values
figure(1)
plot(diag(S),'.','Color',[1 69/255 79/255],'MarkerSize',30)
ylabel('$\sigma_j$','Interpreter','Latex')
title('Singular Values of $\Gamma(X,\dot X)$','Interpreter','Latex')
set(gca,'FontSize',16,'Xlim',[0.9 length(diag(S))+0.1])

% Check which linear combination of nullspace vectors gives known conserved
% quantities:
xi1 = [V(:,end-1) V(:,end)]\[0;0;0;1;0;0;1;0;1; zeros(9,1)]/2; % angular momentum
xi2 = [V(:,end-1) V(:,end)]\[0;0;0;1;0;0;2;0;3; zeros(9,1)]/2; % kinetic energy (Hamiltonian)

%% Plot trajectories on the unit sphere

% Generate unit sphere
[Xs, Ys, Zs] = sphere(100);

% Generate initial conditions on the unit sphere
t = 0:dt:50;
x0 = [-1,-.5,0];
x0 = x0/norm(x0);
[~,xUnit(1,:,:)] = ode45(@(t,x) consODE(x),t,x0,options);
x0 = [-1,0.9,-0.9];
x0 = x0/norm(x0);
[~,xUnit(2,:,:)] = ode45(@(t,x) consODE(x),t,x0,options);
x0 = [0.1,-1,0.3];
x0 = x0/norm(x0);
[~,xUnit(3,:,:)] = ode45(@(t,x) consODE(x),t,x0,options);
x0 = [0.1,-1,-0.3];
x0 = x0/norm(x0);
[~,xUnit(4,:,:)] = ode45(@(t,x) consODE(x),t,x0,options);
x0 = [0.18,-0.9,-0.1];
x0 = x0/norm(x0);
[~,xUnit(5,:,:)] = ode45(@(t,x) consODE(x),t,x0,options);
x0 = [0.3,0.3,0.80];
x0 = x0/norm(x0);
[~,xUnit(6,:,:)] = ode45(@(t,x) consODE(x),t,x0,options);

% Plot trajectories on the unit sphere
figure(2)
plot3(xUnit(1,:,1),xUnit(1,:,2),xUnit(1,:,3),'Color',[1 69/255 79/255],'LineWidth',2)
hold on
plot3(xUnit(2,:,1),xUnit(2,:,2),xUnit(2,:,3),'Color',[146/255 208/255 80/255],'LineWidth',2)
plot3(xUnit(3,:,1),xUnit(3,:,2),xUnit(3,:,3),'Color',[36/255 122/255 254/255],'LineWidth',2)
plot3(xUnit(4,:,1),xUnit(4,:,2),xUnit(4,:,3),'Color',[254/255 174/255 0/255],'LineWidth',2)
plot3(xUnit(5,:,1),xUnit(5,:,2),xUnit(5,:,3),'Color',[255/255 66/255 161/255],'LineWidth',2)
plot3(xUnit(6,:,1),xUnit(6,:,2),xUnit(6,:,3),'Color',[22/255 231/255 207/255],'LineWidth',2)
surf(Xs,Ys,Zs)
colormap('gray')
shading interp
grid off
axis off
axis equal

%% Plot trajectories on the unit ellipsoid

% Generate unit sphere
[Xe, Ye, Ze] = ellipsoid(0,0,0,1,1/sqrt(2),1/sqrt(3),100);

% Generate initial conditions on the unit sphere
t = 0:dt:50;
x0 = [0,-.2,0.1];
x0(1) = -sqrt(1 - 2*x0(2)^2 - 3*x0(3)^2);
[~, xEllip(1,:,:)] = ode45(@(t,x) consODE(x),t,x0,options);
x0 = [0,-.27,0.48];
x0(1) = -sqrt(1 - 2*x0(2)^2 - 3*x0(3)^2);
[~, xEllip(2,:,:)] = ode45(@(t,x) consODE(x),t,x0,options);
x0 = [0,-.63,0.01];
x0(1) = -sqrt(1 - 2*x0(2)^2 - 3*x0(3)^2);
[~, xEllip(3,:,:)] = ode45(@(t,x) consODE(x),t,x0,options);
x0 = [0,-.14,0.55];
x0(1) = -sqrt(1 - 2*x0(2)^2 - 3*x0(3)^2);
[~, xEllip(4,:,:)] = ode45(@(t,x) consODE(x),t,x0,options);
x0 = [0,-.66,-0.19];
x0(1) = -sqrt(1 - 2*x0(2)^2 - 3*x0(3)^2);
[~, xEllip(5,:,:)] = ode45(@(t,x) consODE(x),t,x0,options);
x0 = [0,-.65,0.1];
x0(1) = sqrt(1 - 2*x0(2)^2 - 3*x0(3)^2);
[~, xEllip(6,:,:)] = ode45(@(t,x) consODE(x),t,x0,options);

% Plot trajectories on the unit sphere
figure(3)
plot3(xEllip(1,:,1),xEllip(1,:,2),xEllip(1,:,3),'Color',[1 69/255 79/255],'LineWidth',2)
hold on
plot3(xEllip(2,:,1),xEllip(2,:,2),xEllip(2,:,3),'Color',[146/255 208/255 80/255],'LineWidth',2)
plot3(xEllip(3,:,1),xEllip(3,:,2),xEllip(3,:,3),'Color',[36/255 122/255 254/255],'LineWidth',2)
plot3(xEllip(4,:,1),xEllip(4,:,2),xEllip(4,:,3),'Color',[254/255 174/255 0/255],'LineWidth',2)
plot3(xEllip(5,:,1),xEllip(5,:,2),xEllip(5,:,3),'Color',[255/255 66/255 161/255],'LineWidth',2)
plot3(xEllip(6,:,1),xEllip(6,:,2),xEllip(6,:,3),'Color',[22/255 231/255 207/255],'LineWidth',2)
surf(Xe,Ye,Ze)
colormap('gray')
shading interp
grid off
axis off
axis equal

%% Conservative ODE right-hand-side

function dx = consODE(x)

    dx = [x(2)*x(3); -2*x(1)*x(3); x(1)*x(2)];
    
end