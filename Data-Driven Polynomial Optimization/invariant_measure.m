% -------------------------------------------------------------------------
% Identifying invariant measures from Data
%
% The goal of the script is to use the approximated Lie derivative from
% data to identify invariant measures of the underlying system. We apply
% this method to the chaotic quadratic mapping
%
%           x_{n+1} = 2x_n^2 - 1
%
% and identify both physical and atomic measures. 
% 
% The saved MATLAB file AtomicMeasures.mat contains approximations of the
% atomic measure localized the fixed point x = -0.5 for maxDp = 5,10,20.
% 
% To run this notebook one requires paths to the freely available software
% packages YALMIP and MOSEK. YALMIP can be downloaded at:
%           https://yalmip.github.io/download/
% and MOSEK can be downloaded at:
%           https://www.mosek.com/downloads/
%
% This script accompanies Section 4.5 of Data-Driven Methods for
% Dynamic Systems.
%
% Written by J. Bramburger and G. Fantuzzi.
%
% -------------------------------------------------------------------------

% Clean workspace
clear; 
close all; 
clc
yalmip clear
format long

%% Method Parameters 
% maxDp = max degree of Dp dictionary of obserables (alpha parameter)
% maxDq = max degree of Dq dictionary of obserables (beta parameter)
maxDp = 5;
maxDq = 2*maxDp;
N = 1e4;
TOL = 0;

%% Generate synthetic data

% Chebfun objects
T1 = chebpoly(0:maxDp,[-1,1]);
T2 = chebpoly(0:maxDq,[-1,1]);

% Generate data
x = zeros(N,1);
x(1) = 0.25;
for n = 2:N
    % iterating the map 2x^2 - 1
    x(n) = 2*x(n-1)^2 - 1; 
end

% EDMD matrix
Q = T2( x(1:N-1) )'; 
P = T1( x(2:N) )'; 
K = edmd_with_thresholding(Q,P,TOL);

% Lie Derivative
L = K - eye(size(K));

%% Moment Matrices

% Moment vector
y = sdpvar(maxDq,1);

% Moment matrix
A = chebsdp_1d(maxDp);
M0 = reshape(A*[1;y],[maxDp+1,maxDp+1]);

% Localizing moment matrix for [-1,1]
B = chebsdp_1d_locball(maxDp);
M1 = reshape(B*[1;y],[maxDp,maxDp]);

%% Optimization procedure to identify physical meausres

% Approximate first moment from data
yObj = sum(x)/N;
%yObj = 0; % <--- exact moment value

% Solve
OBJ = ( y(1) - yObj ).^2;
sol = optimize([L*[1;y]==0, M0>=0, M1>=0], OBJ)

%% Build approximation to density with chebfun

M = A(:,1) * 2;
for j = 1:maxDq
    T = chebpoly(j,[-1,1]);
    M = M + A(:,j+1) .* sum(T);
end

M = full(reshape(M,[maxDp+1,maxDp+1]));
rho = chebfun(M\value([1;y(1:maxDp)]),'coeffs');
xx = linspace(-1,1,10000);

% Plot results
figure(1)
plot(xx, 1./pi./sqrt(1-xx.^2),'--','Color',[1 69/255 79/255],'LineWidth',3)
hold on
plot(rho,'Color',[36/255 122/255 254/255],'LineWidth',4)
axis([-1 1 0 4])
xlabel('$x$','interpreter','latex')
legend('Exact','Discovered','Location','Best','interpreter','latex','FontSize',20)
set(gca,'FontSize',16)

%% Optimization procedure for discovering (approximately) atomic measures

% Solve with linear objective function
OBJ = y(1);
sol = optimize([L*[1;y]==0, M0>=0, M1>=0], OBJ)

%% Build approximation to density with chebfun and plot it

% Get density
M = A(:,1) * 2;
for j = 1:maxDq
    T = chebpoly(j,[-1,1]);
    M = M + A(:,j+1) .* sum(T);
end
M = full(reshape(M,[maxDp+1,maxDp+1]));
rho = chebfun(M\value([1;y(1:maxDp)]),'coeffs');

% Space points
xx = linspace(-1,1,10000);

% Plot result
figure(1)
plot(xx,rho(xx),'Color',[1 69/255 79/255],'LineWidth',2);

%% Print minizers
% --> uses cheb_extractMinimizers.m which is a modified extractMinimizers 
%      for the Chebyshev basis 

xopt = cheb_extractMinimizers(value(M0), 0:maxDq/2)
hold on
plot(xopt,rho(xopt),'ko','MarkerSize',8,'MarkerFaceColor','k')

%% Plot atomic measures results

load AtomicMeasures.mat

figure(3)
plot(xx,measPltal5,'Color',[0 168/255 0],'LineWidth',4)
hold on 
plot(xx,measPltal10,'Color',[1 69/255 79/255],'LineWidth',4)
plot(xx,measPltal20,'Color',[36/255 122/255 254/255],'LineWidth',4)
xlabel('$x$','interpreter','latex')
legend('$\alpha = 5$','$\alpha = 10$','$\alpha = 20$','Location','Best','interpreter','latex','FontSize',20)
set(gca,'FontSize',16)

