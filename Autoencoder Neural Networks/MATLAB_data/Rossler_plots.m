% -------------------------------------------------------------------------
% The purpose of this script is to make the plots corresponding to the
% trained autoencoder neural networks for the Rossler Poincare map. The 
% plots take in data from the jupyter notebook Rossler_conj.ipynb and plot
% the results for presentation. 
%
% This script accompanies Section 6.5 of Data-Driven Methods for
% Dynamic Systems. 
%
% Author: Jason J. Bramburger
% -------------------------------------------------------------------------

%% Poincare Section and Mapping Plots 

% Clean workspace
clear all; close all; clc

%Model parameters 
a = 0.1;
b = 0.1;
mu = 18; 

%ODE generation parameters
m = 3; %Dimension of ODE
n = m-1; %Dimension of Poincare section
dt = 0.001;
tspan = (0:10000000-1)*dt;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,m));

x0(1,:) = [0.0; -15; 0]; 
[~,xdat]=ode45(@(t,x) Rossler(x,a,b,mu),tspan,x0(1,:),options);

%Initialize
Psec = [];
PsecNext = [];

%Create Poincare section data
temp = [];
count = 1;
for j = 1:length(xdat(:,1))-1 
    if  (xdat(j,1) < 0.1 && xdat(j+1,1) >= 0.1) %&& j >= 0.5*length(xdat(i,:,1))) 
        temp(count) = xdat(j+1,2); %nth iterate
        count = count + 1;
    end
end
Psec = temp(1:length(temp)-1);
PsecNext = temp(2:length(temp));

%% Plot Poincare section data 

figure(1)
plot(Psec,PsecNext,'k.','MarkerSize',20)
xlabel('$x_2$','Interpreter','Latex')
ylabel('$f(x_2)$','Interpreter','Latex')
%set(gca,'FontSize',16,'Xlim',[-14.2 -9.7],'Ylim',[-14.2 -9.7]) % for mu = 9
set(gca,'FontSize',16,'Xlim',[-27 -11.5],'Ylim',[-27 -11.5]) % for mu = 11, 13, 18
box on


%% Comparing c = 18 results

% Clean workspace
clear all; close all; clc

load Rossler_data.mat

% Discovered SINDy map
slope = 0.06734752269044327; % used for scaling data
yint = 1.791773357461594;
xsindy = 0.9776 - 3.5219*(slope*xn + yint) + 3.4089*(slope*xn + yint).^2;
xsindy = (xsindy - yint)/slope;

% Plot Poincare map reconstruction with quadratic fit
figure(1)
hold on
plot(xn,xnp1,'k.','MarkerSize',10)
plot(xn,xsindy,'.','Color',[36/255 122/255 254/255],'MarkerSize',10)
xlabel('$x_2$','Interpreter','Latex')
ylabel('$f(x_2)$','Interpreter','Latex')
legend('Section Data','Quadratic Fit','Interpreter','Latex','Location','Best','FontSize',20)
set(gca,'FontSize',16,'Xlim',[-27 -12],'Ylim',[-27 -12])
box on

% Plot Poincare map reconstruction with autoencoder
figure(2)
hold on
plot(xn,xnp1,'k.','MarkerSize',10)
plot(xn,xnp1_recon,'.','Color',[1 69/255 79/255],'MarkerSize',10)
xlabel('$x_2$','Interpreter','Latex')
ylabel('$f(x_2)$','Interpreter','Latex')
legend('Section Data','Autoencoder','Interpreter','Latex','Location','Best','FontSize',20)
set(gca,'FontSize',16,'Xlim',[-27 -12],'Ylim',[-27 -12])
box on

% Plot forcasts between autoencoder and quadratic map

% x-axis data
n = 0:length(xpred)-1;

% SINDy step
xsindypred = zeros(length(n),1);
xsindypred(1) = xn(1);
for ind = 2:length(n)
    xsindypred(ind) = 0.9776 - 3.5219*(slope*xsindypred(ind-1) + yint) + 3.4089*(slope*xsindypred(ind-1) + yint)^2; 
    xsindypred(ind) = (xsindypred(ind) - yint)/slope;
end

xpred = [xn(1); xpred];

figure(3)
clf
hold on
plot(n,xn(1:n(end)+1),'k.-','MarkerSize',20,'LineWidth',2)
plot(n,xsindypred,'.-','Color',[36/255 122/255 254/255],'MarkerSize',20,'LineWidth',2)
plot(n,xpred(1:n(end)+1),'.-','Color',[1 69/255 79/255],'MarkerSize',20,'LineWidth',2)
xlabel('Iterates','Interpreter','Latex')
ylabel('Forecast','Interpreter','Latex')
set(gca,'FontSize',16,'Xlim',[0,10])
box on

%% Rossler right-hand-side

function dx = Rossler(x,a,b,c)

    dx = [-x(2) - x(3); x(1) + a*x(2); b + x(3)*(x(1) - c)];

end








