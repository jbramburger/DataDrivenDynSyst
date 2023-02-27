% -------------------------------------------------------------------------
% Plotting upper and lower bounds on homoclinic existence
%
% This script plots the upper and lower bounds over a range of the
% hyperparameter lambda for the existence of homoclinic orbit. Data is
% loaded in from the .mat files homoclinic_lambda_ + 'upper' or 'lower' for
% upper or lower bounds + 'm2' or 'm3' for m = 2 or 3, respectively. Each
% file containts the lambda values, the value of m for reference, and
% arrays c2, c3, c4, c5 representing the values of upper or lower bounds 
% on mu for degree 2,3,4, or 5 auxiliary functions, respectively. 
%
% This script accompanies Section 4.2 of Data-Driven Methods for
% Dynamic Systems. 
%
% Author: Jason J. Bramburger
% -------------------------------------------------------------------------

% Clean workspace
clear all; 
close all;
clc

%% Load in upper bound data

% m = 2 data
load homoclinic_lambda_upper_m2.mat
load homoclinic_lambda_upper_m2_part2.mat

% m = 3 data
load homoclinic_lambda_upper_m3.mat
load homoclinic_lambda_upper_m3_part2.mat

%% Unpack data and aggregate

lambda = [lambda lambda_2];
c3 = [c3 c3_2];
c4 = [c4 c4_2];
c5 = [c5(1:181) c5_2];

% Move failures off the figure
if m == 2 
    c2 = [c2 2*ones(1,length(lambda_2))];
    c2(1:46) = 2;
    c2(174:end) = 2;
    c3(260:end) = 2;
else
   c2 = [c2 c2_2];
end


%% Plot Upper bounds

figure(1)
hold on
plot(lambda(1:3:end),c2(1:3:end),'Color',[1 69/255 79/255],'LineWidth',4)
hold on
plot(lambda(1:3:end),c3(1:3:end),'--','Color',[36/255 122/255 254/255],'LineWidth',4)
plot(lambda(1:3:end),c4(1:3:end),'-.','Color',[0 168/255 0],'LineWidth',4)
plot(lambda(1:3:end),c5(1:3:end),':','Color',[255/255 66/255 161/255],'LineWidth',4)

%% Load lower bound data
if m == 2
    load homoclinic_lambda_lower_m2.mat
elseif m == 3
    load homoclinic_lambda_lower_m3.mat
end

%% Plot lower bounds

if m == 2 % m = 2 Specifications
    figure(1)
    hold on
    plot(lambda(7:37),c3(7:37),'--','Color',[36/255 122/255 254/255],'LineWidth',4)
    plot(lambda(7:37),c4(7:37),'-.','Color',[0 168/255 0],'LineWidth',4)
    plot(lambda(7:2:39),c5(7:2:39),':','Color',[255/255 66/255 161/255],'LineWidth',4)
    set(gca,'FontSize',14)
    axis([0 lambda(37) 0.5 1])
    yticks(0.5:0.1:1)
elseif m == 3 % m = 3 Specifications
    figure(1)
    hold on
    plot(lambda(7:end),c3(7:end),'--','Color',[36/255 122/255 254/255],'LineWidth',4)
    plot(lambda(7:end),c4(7:end),'-.','Color',[0 168/255 0],'LineWidth',4)
    plot(lambda(7:2:end),c5(7:2:end),':','Color',[255/255 66/255 161/255],'LineWidth',4)
    set(gca,'FontSize',14)
    axis([0 lambda(37) 0.30 0.7])
    yticks(0.3:0.1:0.765)
end

%legend({'d = 2','d = 3', 'd = 4','d = 5'}, 'Interpreter','latex','Location','best')
set(gca,'FontSize',20)
xlabel('$\lambda$','Interpreter','latex','FontSize',24,'FontWeight','Bold')
ylabel('$\mu$','Interpreter','latex','FontSize',24,'FontWeight','Bold')
box on








