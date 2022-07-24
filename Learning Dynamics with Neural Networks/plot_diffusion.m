% -------------------------------------------------------------------------
% Plotting the PINN Solution to the Diffusion Equation
%
% This script accompanies the jupyter notebook Diffusion_PINN.ipynb. Here
% we take the PINN solution to the diffusion equation from the notebook,
% plot it and compare the final state at t = 1 to the true final state
% exp(-0.1*pi^2)*cos(pi*x).
%
% To recreate the results from the text simply load in diffusion_sol.mat. 
%
% This script accompanies Section 4.4 of Data-Driven Methods for
% Dynamic Systems. 
%
% Author: Jason J. Bramburger
% -------------------------------------------------------------------------

% Clean workspace
clear all; close all; clc

%% Load in PINN data
load diffusion_sol.mat

%% Plot PINN solution

% For plotting initial and final solutions
x = X(:,1);
n = length(x);

% Plot PINN solution
figure(1)
hold on
surf(T,X,U)
plot3(0*x,x,U(:,1),'Color',[1 69/255 79/255],'LineWidth',5)
plot3(ones(1,n),x,U(:,end),'Color',[36/255 122/255 254/255],'LineWidth',5)
hold off
grid on
view(-45,25)
shading interp
colormap gray
xlabel('t')
ylabel('x')
zlabel('u(x,t)')
set(gca,'FontSize',16,'Xlim',[0 1],'Ylim',[-1 1])

%% Compare final state to true final state

ufinal = exp(-0.1*pi^2)*cos(pi*x);

figure(2)
plot(x,U(:,end),'Color',[36/255 122/255 254/255],'LineWidth',5)
hold on
plot(x,ufinal,'k--','LineWidth',3)
set(gca,'FontSize',16,'Xlim',[-1 1])
legend('PINN Solution','True Solution')

% Compute sup error
error = norm(U(:,end) - ufinal,'inf')




