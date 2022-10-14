% -------------------------------------------------------------------------
% The purpose of this script is to make the plots corresponding to the
% trained autoencoder neural networks. The plots take in data from the
% corresponding jupyter notebooks and plot the results for presentation.
%
% To recreate the results from the text simply load in diffusion_sol.mat. 
%
% This script accompanies Section 6 of Data-Driven Methods for
% Dynamic Systems. 
%
% Author: Jason J. Bramburger
% -------------------------------------------------------------------------

% Clean workspace
clear all; close all; clc

%% Tent --> Logistic Conjugacy results 
%   Date comes from Tent2Logistic.ipynb notebook

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

%% Tent --> Sine Conjugacy results 
%   Date comes from Tent2Sine.ipynb notebook

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















