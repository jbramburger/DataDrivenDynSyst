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



