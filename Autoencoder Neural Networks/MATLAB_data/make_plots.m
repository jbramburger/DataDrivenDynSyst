% -------------------------------------------------------------------------
% The purpose of this script is to make the plots corresponding to the
% trained autoencoder neural networks. The plots take in data from the
% corresponding jupyter notebooks and plot the results for presentation. 
%
% This script accompanies Chapter 6 of Data-Driven Methods for
% Dynamic Systems. 
%
% Author: Jason J. Bramburger
% -------------------------------------------------------------------------


%% Tent --> Logistic Conjugacy results 
%   Data comes from Tent2Logistic.ipynb notebook

% Clean workspace
clear all; close all; clc

load Tent2Logistic_Conj.mat 
% --> Outputs are (x,y,xdecode)
%        x = linspace(0,1,1000)
%        y = encoder(x) - approximate conjugacy output
%        xdecode = decoder(x) - approximate inverse conjugacy output

% Plot encoder & conjugacy to compare
h = sin(pi*x/2).^2; % exact conjugacy

figure(1)
hold on
plot(x,y,'Color',[36/255 122/255 254/255],'LineWidth',6)
plot(x,h,'k--','LineWidth',3)
xlabel('$x$','Interpreter','Latex')
ylabel('$y$','Interpreter','Latex')
legend('$\eta(x)$','$h(x)$','Interpreter','Latex','Location','Best','FontSize',20)
set(gca,'FontSize',16,'Xlim',[0 1],'Ylim',[0 1])
box on

% Plot decoder & inverse conjugacy
hinv = 2*asin(sqrt(x))/pi; % exact inverse of conjugacy

figure(2)
hold on
plot(x,xdecode,'Color',[1 69/255 79/255],'LineWidth',6)
plot(x,hinv,'k--','LineWidth',3)
xlabel('$y$','Interpreter','Latex')
ylabel('$x$','Interpreter','Latex')
legend('$\rho(x)$','$h^{-1}(x)$','Interpreter','Latex','Location','Best','FontSize',20)
set(gca,'FontSize',16,'Xlim',[0 1],'Ylim',[0 1])
box on

% Pointwise errors in plotted points
conjError = max(abs(y - h'))
invconjError = max(abs(xdecode - hinv'))

%% Tent --> Sine Conjugacy results 
%   Data comes from Tent2Sine.ipynb notebook

% Clean workspace
clear all; close all; clc

load Tent2Sine_Conj.mat 
% --> Outputs are (x,y,xdecode)
%        x = linspace(0,1,1000)
%        y = encoder(x) - approximate conjugacy output
%        xdecode = decoder(x) - approximate inverse conjugacy output

% Plot encoder 
figure(1)
hold on
plot(x,y,'Color',[36/255 122/255 254/255],'LineWidth',6)
xlabel('$x$','Interpreter','Latex')
ylabel('$y$','Interpreter','Latex')
set(gca,'FontSize',16,'Xlim',[0 1],'Ylim',[-1 0])
box on

% Plot decoder
figure(2)
hold on
plot(-x,xdecode,'Color',[1 69/255 79/255],'LineWidth',6)
xlabel('$y$','Interpreter','Latex')
ylabel('$x$','Interpreter','Latex')
set(gca,'FontSize',16,'Xlim',[-1 0],'Ylim',[0 1])
box on

%% Saddle-node normal form branch reconstruction

load saddle_node.mat
% --> Outputs are (fixed_s,fixed_u)
%        stable (s) and unstable (u) branches of bifurcating solutions

% Taylor expansion approximation
mu = 0:0.0001:max(fixed_s(:,2));
xTaylor_s = sqrt(2*mu/exp(1));
xTaylor_u = -sqrt(2*mu/exp(1));

% Find exact fixed points
for j = 1:length(mu)
   snEqn = @(x) double(mu(j)+exp(1))*x - exp(x); 
   xp(j) = fsolve(snEqn, 1.5);
   xm(j) = fsolve(snEqn, 0.5);
end

% Bifurcation diagram with Taylor approximation on it
figure(1)
plot(mu + exp(1),xTaylor_s + 1,'k','LineWidth',2)
hold on
plot(mu + exp(1),xTaylor_u + 1,'k--','LineWidth',2)
plot(fixed_s(:,2) + exp(1),fixed_s(:,1) + 1,'Color',[1 69/255 79/255],'LineWidth',3)
plot(fixed_u(:,2) + exp(1),fixed_u(:,1) + 1,'--','Color',[1 69/255 79/255],'LineWidth',3)
plot(exp(1),1,'.','Color',[36/255 122/255 254/255],'MarkerSize',40)
xlabel('$\mu$','Interpreter','Latex')
ylabel('$x^*(\mu)$','Interpreter','Latex')
set(gca,'FontSize',16,'Xlim',[exp(1) - 0.02 3.01])
box on

% (Numerically) exact bifurcation diagram
figure(2)
plot(fixed_s(:,2) + exp(1),fixed_s(:,1) + 1,'Color',[1 69/255 79/255],'LineWidth',5)
hold on
plot(mu + exp(1),xp,'k--','LineWidth',2)
plot(fixed_u(:,2) + exp(1),fixed_u(:,1) + 1,'Color',[1 69/255 79/255],'LineWidth',5)
plot(mu + exp(1),xm,'k--','LineWidth',2)
plot(exp(1),1,'.','Color',[36/255 122/255 254/255],'MarkerSize',40)
xlabel('$\mu$','Interpreter','Latex')
ylabel('$x^*(\mu)$','Interpreter','Latex')
set(gca,'FontSize',16,'Xlim',[exp(1) - 0.02 3.01])
box on

%% Period-doubling normal form bifurcation diagram

% Clean workspace
close all
clear all
clc

% Initialization for bifurcation diagram
N = 50000; % number of iterations
mu = 0.22:0.0001:0.35;
x = zeros(N,3);
bifDiag = zeros(length(mu),10);

% Generate data for bifurcation diagram
for j = 1:length(mu)
    x(1,:) = [1; 1; -mu(j)] + 0.1;
    for n = 1:N-1
        
        x(n+1,1) = x(n,2);
        x(n+1,2) = 1 + mu(j)*x(n,1) + x(n,2)*x(n,3);
        x(n+1,3) = x(n,3) - mu(j)*x(n,1)*x(n,2) - x(n,3)*x(n,2)^2;
        
    end
    
    bifDiag1(j,:) = x(end-9:end,1) - 1;
    bifDiag3(j,:) = x(end-9:end,3) + mu(j);
end

% Plot attractors in a bifurcation diagram
figure(1)
plot(mu - 0.25,bifDiag1,'k.','MarkerSize',10)
xlabel('$\mu - 0.25$','Interpreter','Latex')
ylabel('$x_3 - \mu$','Interpreter','Latex')
set(gca,'FontSize',16,'Xlim',[mu(1)-0.25, mu(end)-0.25],'Ylim',[min(bifDiag1(end,:)), max(bifDiag1(end,:))])
box on

% Load normal form data
load period_doubling.mat
% --> Outputs are (flip_x1,flip_x2,fixed_pt_x)
%        2-cycle branches of bifurcating solutions are given in flip_x1 and
%        flip_x2, while the trivial/fixed point branch is given in
%        fixed_pt_x. All outputs are 4-dimensional: 3 state variables and
%        one parameter.

% Plot bifurcation diagram from normal form
figure(2)
plot(fixed_pt_x(1:10000,4)/10,fixed_pt_x(1:10000,1),'k','LineWidth',3)
hold on
plot(fixed_pt_x(10001:20000,4)/10,fixed_pt_x(10001:20000,1),'k--','LineWidth',3)
plot(flip_x1(:,4)/10,flip_x1(:,1),'k','LineWidth',3)
plot(flip_x2(:,4)/10,flip_x2(:,1),'k','LineWidth',3)
plot(mu - 0.25,bifDiag1,'b.','MarkerSize',10)
xlabel('$\mu$','Interpreter','Latex')
ylabel('$x_1$','Interpreter','Latex')
%set(gca,'FontSize',16,'Xlim',[0.23, 0.27])
box on

% x3 bifurcation diagram 
figure(3)
plot(fixed_pt_x(1:10000,4)/10 + 0.25,fixed_pt_x(1:10000,3) - fixed_pt_x(1:10000,4)/10 - 0.25,'k','LineWidth',3)
hold on
plot(fixed_pt_x(10001:20000,4)/10 + 0.25,fixed_pt_x(10001:20000,3) - fixed_pt_x(10001:20000,4)/10 - 0.25,'k--','LineWidth',2)
plot(flip_x1(:,4)/10 + 0.25,flip_x1(:,3) - flip_x1(:,4)/10 - 0.25,'k','LineWidth',3)
plot(flip_x2(:,4)/10 + 0.25,flip_x2(:,3) - flip_x2(:,4)/10 - 0.25,'k','LineWidth',3)
%plot(mu - 0.25,bifDiag3,'b.','MarkerSize',10)
xlabel('$\mu$','Interpreter','Latex')
ylabel('$x_3$','Interpreter','Latex')
set(gca,'FontSize',16,'Xlim',[0.23, 0.27])
box on

%% Koopman eigenfunctions for global linearization 
%   Data comes from Global_Linearization.ipynb notebook

% Clean workspace
clear all; close all; clc

%load koopObs_fixed.mat     %   <---- fixed/static eigenvalues
load koopObs_variable.mat %   <---- variable eigenvalues
% --> Outputs are (x,koop)
%        x = 16000 evenly spaced grid points in [-2,2]x[-2,2]
%        y = encoder(x) - value of each of the 3 Koopman observables at
%        each x

[X1,X2] = meshgrid(-2:0.01:2,-2:0.01:2);

% Plot Discovered Koopman Eigenfunctions
figure(1)
set(gcf, 'Position',  [300, 400, 700, 200])

% Eigenfunction 1
subplot(1,3,1)
surf(X1,X2,reshape(koop(:,1),[401 401]))
shading interp
view(0,90)
set(gca,'FontSize',14,'Xlim',[-2,2],'Ylim',[-2,2])

% Eigenfunction 2
subplot(1,3,2)
surf(X1,X2,reshape(koop(:,2),[401 401]))
shading interp
view(0,90)
set(gca,'FontSize',14,'Xlim',[-2,2],'Ylim',[-2,2])

% Eigenfunction 3
subplot(1,3,3)
surf(X1,X2,reshape(koop(:,3),[401 401]))
shading interp
view(0,90)
set(gca,'FontSize',14,'Xlim',[-2,2],'Ylim',[-2,2])

%% Koopman eigenfunctions for global linearization (Complex eigenvalue case) 
%   Data comes from Global_Linearization.ipynb notebook

% Clean workspace
clear all; close all; clc

%load koopObs_fixed.mat     %   <---- fixed/static eigenvalues
load koopObs_complex.mat %   <---- variable eigenvalues
% --> Outputs are (x,koop)
%        x = 16000 evenly spaced grid points in [-2,2]x[-2,2]
%        y = encoder(x) - value of each of the 6 (3 real + 3 imaginary part) 
%           Koopman observables at each x

[X1,X2] = meshgrid(-2:0.01:2,-2:0.01:2);

% Plot Discovered Koopman Eigenfunctions
figure(1)
set(gcf, 'Position',  [300, 400, 700, 200])

% Exact Eigenfunction 1
subplot(1,3,1)
surf(X1,X2,reshape(X1,[401 401]))
shading interp
view(0,90)
set(gca,'FontSize',14,'Xlim',[-2,2],'Ylim',[-2,2])

% Eigenfunction 1 "Real Part"
subplot(1,3,2)
surf(X1,X2,reshape(koop(:,1),[401 401]))
shading interp
view(0,90)
set(gca,'FontSize',14,'Xlim',[-2,2],'Ylim',[-2,2])

% Eigenfunction 1 "Imaginary Part"
subplot(1,3,3)
surf(X1,X2,reshape(koop(:,2),[401 401]))
shading interp
view(0,90)
set(gca,'FontSize',14,'Xlim',[-2,2],'Ylim',[-2,2])

%% Torus plot

%Create R and THETA data
theta = 0:pi/100:2*pi;
r = 2*pi:pi/200:3*pi;
[R,T] = meshgrid(r,theta);

%Create top and bottom halves
Z_top = 2*sin(R);
Z_bottom = -2*sin(R);

% Add in winding torus trajectory
t = 80:0.1:300;
x = (2.5*pi + 0.5*pi*cos(t)).*cos(t/pi);
y = (2.5*pi + 0.5*pi*cos(t)).*sin(t/pi);
z = 0.5*pi*sin(t);
    
%Convert to Cartesian coordinates and plot
[X,Y,Z] = pol2cart(T,R,Z_top);
surf(X,Y,Z);
hold on;
[X,Y,Z] = pol2cart(T,R,Z_bottom);
surf(X,Y,Z);
colormap gray
alpha 0.8
shading interp
axis equal
axis off

% Plot trajectory on torus
k = 230;
plot3(x(1:k),y(1:k),z(1:k),'Color',[1 69/255 79/255],'LineWidth',2)

%% Plot Action-Angle Omega Neural Network Output
%   Data comes from ActionAngle.ipynb notebook

% Clean workspace
clear all; close all; clc
  
load Omega_Network.mat    
% --> Outputs are (r,omega)
%        r = 1000 evenly spaced points in [0,1]
%        omega = Omega(r) - value of the trained neural network that
%                identifies the Koopman eigenvalue omega(r)

% Plot Omega neural network output
figure(1)
plot(r,omega,'Color',[0/255 168/255 0/255],'LineWidth',5)
hold on
plot(r(499),omega(499),'.','Color',[1 69/255 79/255],'MarkerSize',50)
xlabel('$y_1^2 + y_2^2$','Interpreter','Latex')
ylabel('$\Omega(y_1^2 + y_2^2)$','Interpreter','Latex')
set(gca,'FontSize',16,'Xlim',[0 0.3],'Ylim',[0.60 1.01])
box on

% Generate sample data from the Kepler system

M = 7; % number of trajectories
x0 = [linspace(1,2,M); zeros(1,M)]';
dt = 0.02;
t = 0:dt:10;
x = zeros(length(t),2,M);

% Plot trajectories of the Kepler problem
figure(2)
hold on
for m = 1:M 
    sol = [];
    [t, sol] = ode45(@(t,x) Kepler(t,x),0:dt:(10+4/m),x0(m,:));
    if m == (M+1)/2
        plot(sol(:,1),sol(:,2),'Color',[1 69/255 79/255],'LineWidth',3)
    else
        plot(sol(:,1),sol(:,2),'Color',[0.5 0.5 0.5],'LineWidth',3)
    end
end
plot(1,0,'.','Color',[36/255 122/255 254/255],'MarkerSize',30)
set(gca,'FontSize',16,'Xlim',[0.6 2.05],'Ylim',[-0.55 0.55])
box on

% Plot trajectories of the linearized system
th = 0:0.01:2*pi;
rad = linspace(0,0.3,M);

figure(3)
hold on
for m = 1:M
   sol = [];
   sol = sqrt(rad(m))*[cos(th); sin(th)]';
   if m == (M+1)/2
        plot(sol(:,1),sol(:,2),'Color',[1 69/255 79/255],'LineWidth',3)
    else
        plot(sol(:,1),sol(:,2),'Color',[0.5 0.5 0.5],'LineWidth',3)
    end
end
plot(0,0,'.','Color',[36/255 122/255 254/255],'MarkerSize',30)
set(gca,'FontSize',16,'Xlim',[-sqrt(0.35) sqrt(0.35)],'Ylim',[-sqrt(0.35) sqrt(0.35)])
box on

%% Kepler RHS

function rhs = Kepler(t,x)
    rhs = [x(2); x(1)^(-3) - x(1)^(-2)];
end


