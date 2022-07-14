% -------------------------------------------------------------------------
% Plotting the Henon Neural Network Prediction 
%
% This script accompanies the jupyter notebook Forecast.ipynb. Here
% we take the Henon prediction from the feedforward neural network and 
% compare it to the true iterates of the Henon map.
%
% A number of different results can be loaded using the main variation
%       Henon_step=#_b=0_&.mat
% where # can be 1,2,3, or 4 and & can be 0.01 or 0.3.
%
% This script accompanies Section 4.3 of Data-Driven Methods for
% Dynamic Systems. 
%
% Author: Jason J. Bramburger
% -------------------------------------------------------------------------

% Clean workspace
clear all; close all; clc


%% Load Henon Prediction data
load Henon_step=4_b=0_3

%% Plot Prediction vs Truth

% Compare x_1 predictions
figure(1)
clf
plot(xtrue(1:30,1),'k.-','MarkerSize',20)
hold on
plot(xpred(1:30,1),'.-','Color',[1 69/255 79/255],'MarkerSize',20)
set(gca,'FontSize',16,'Xlim',[1,30])

% Compare x_2 predictioins
figure(2)
clf
plot(xtrue(1:30,2),'k.-','MarkerSize',20)
hold on
plot(xpred(1:30,2),'.-','Color',[36/255 122/255 254/255],'MarkerSize',20)
set(gca,'FontSize',16,'Xlim',[1,30])


