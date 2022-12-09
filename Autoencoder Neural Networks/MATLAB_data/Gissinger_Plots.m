% -------------------------------------------------------------------------
% The purpose of this script is to make the plots corresponding to the
% trained autoencoder neural networks for the Gissinger Poincare map. The 
% plots take in data from the jupyter notebook Gissinger_conj.ipynb and 
% plot the results for presentation. 
%
% This script accompanies Section 6.5 of Data-Driven Methods for
% Dynamic Systems. 
%
% Author: Jason J. Bramburger
% -------------------------------------------------------------------------

%% Poincare Section and Mapping Plots 

% Clean workspace
clear all; close all; clc

% Gissinger parameters
nu = 0.1;
gam = 0.85;
mu = 0.12;

%ODE generation parameters
m = 3; %Dimension of ODE
n = m-1; %Dimension of Poincare section
dt = 0.001;
tspan = (0:10000000-1)*dt;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,m));

x0(1,:) = [-1.01; 1.01; 1.01]; 
[~,xdat]=ode45(@(t,x) Gissinger(x,mu,nu,gam),tspan,x0(1,:),options);

%Initialize
Psec = [];
PsecNext = [];

%Create Poincare section data
temp = [];
count = 1;
for j = 1:length(xdat(:,1))-1 
    if ((xdat(j,1) + xdat(j,2)) < 0) && ((xdat(j+1,1) + xdat(j+1,2)) >= 0)
        temp(count,:) = xdat(j+1,2:3); %nth iterate
        count = count + 1;
    end
end
Psec = temp(1:length(temp)-1,:);
PsecNext = temp(2:length(temp),:);

%% Plot Gissinger Attractor

figure(1)
plot3(xdat(1:1000000,1),xdat(1:1000000,2),xdat(1:1000000,3),'Color',[36/255 122/255 254/255],'LineWidth',1)
hold on
plot3(-1,1,1,'.','Color',[1 69/255 79/255],'MarkerSize',50)
xlabel('$x_1(t)$','Interpreter','Latex')
ylabel('$x_2(t)$','Interpreter','Latex')
zlabel('$x_3(t)$','Interpreter','Latex')
set(gca,'FontSize',16)
axis tight
box on
grid on


%% Plot Poincare section data 

figure(2)
plot(Psec(:,1),PsecNext(:,1),'.','Color',[0 168/255 0],'MarkerSize',20)
xlabel('$\hat x_2$','Interpreter','Latex')
ylabel('$f_1(\hat x_2,\hat x_3)$','Interpreter','Latex')
set(gca,'FontSize',16)
box on
axis tight

figure(3)
plot(Psec(:,2),PsecNext(:,2),'.','Color',[0 168/255 0],'MarkerSize',20)
xlabel('$\hat x_3$','Interpreter','Latex')
ylabel('$f_2(\hat x_2,\hat x_3)$','Interpreter','Latex')
set(gca,'FontSize',16)
box on
axis tight

%% Plotting Gissinger map

y = linspace(0,1.18,1000);

figure(4)
plot(y,8.53442*y - 18.29989*y.^2 + 9.81721*y.^3,'k','LineWidth',2)
hold on
plot(0,0,'.','Color',[1 69/255 79/255],'MarkerSize',50)
xlabel('$y$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$g(y)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
set(gca,'fontsize',16)
box on
axis tight
grid on


%% Plotting UPOs

% Poincare section line
pline = linspace(1,2,100);

% Generate UPO using Gissinger_UPO.m function
[x1UPO, x2UPO, x3UPO, period] = Gissinger_UPO('1211');
close all % suppress output of Gissinger_UPO

figure(5)
% Equilibrium point
plot(-1,1,'.','Color',[1 69/255 79/255],'MarkerSize',50)
hold on
% attractor in grey
plot(xdat(1:1000000,1),xdat(1:1000000,2),'Color',[0.5 0.5 0.5],'LineWidth',1)
% Line of Poincare section in (x_1,x_2) plane
plot(-pline,pline,'Color',[0 168/255 0],'LineWidth',2)
% Plot UPO
plot([x1UPO; x1UPO],[x2UPO; x2UPO],'Color',[36/255 122/255 254/255],'LineWidth',5)
xlabel('$\hat x(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$\hat y(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
set(gca,'fontsize',16)
box on
axis tight
grid on



%% Gissinger RHS

function dx = Gissinger(x,mu,nu,gamma)
    
    % Equilibria for scaling
    xstar = sqrt(nu + gamma*sqrt(nu/mu));
    ystar = sqrt(mu + gamma*sqrt(mu/nu));
    zstar = -sqrt(nu*mu) - gamma;
    
    % Rescaled variables
    x1hat = x(1)*xstar;
    x2hat = x(2)*ystar;
    x3hat = x(3)*zstar;

    dx = [(mu*x1hat - x2hat*(x3hat + gamma))/xstar; (-nu*x2hat + x1hat*(x3hat + gamma))/ystar; (-x3hat + x1hat*x2hat)/zstar];
    
end