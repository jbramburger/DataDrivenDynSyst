% -------------------------------------------------------------------------
% Identifying Lyapunov Functions from Data with EDMD
%
% The goal of this script is to use EDMD and semidefinite programming to 
% identify a Lyapunov function for the planar discrete-time system
%      x1 --> 0.3*x1 , 
%      x2 --> -x1 + 0.5*x2 + (7/18)*x1^2,
% using only data gathered from the system.
% 
% To run this notebook one requires paths to the freely available software
% packages YALMIP and MOSEK. YALMIP can be downloaded at:
%           https://yalmip.github.io/download/
% and MOSEK can be downloaded at:
%           https://www.mosek.com/downloads/
%
% This script accompanies Section 4.3 of Data-Driven Methods for
% Dynamic Systems.
%
% Written by J. Bramburger and G. Fantuzzi.
%
% -------------------------------------------------------------------------

% Clean workspace
clear; close all; clc
yalmip clear
format long

%% Method Parameters 
% maxDp = max degree of Dp dictionary of obserables
% maxDq = max degree of Dq dictionary of obserables
% epsilon = hyperparameter specific to Lyapunov function for sharp bounds
% cleanVal = remove coefficients smaller than this value from Lyap function
maxDp = 4;
maxDq = maxDp + 2;
epsilon = 1;
cleanVal = 1e-4;

%% Generate synthetic data

% Number of data points
N = 1e2; 

% Random points in [-2,2]x[-2,2]
xdat = 2*rand(N,2) - 1;
ydat = zeros(N,2);

% Images of random points under map dynamics
ydat(:,1) = 0.3*xdat(:,1);
ydat(:,2) = -xdat(:,1) + 0.5*xdat(:,2) + (7/18)*xdat(:,1).^2;

%% Create P and Q matrices

% Q matrix
pow = monpowers(2,maxDq);
ell = size(pow,1); % number of nontrivial monomials in Q
Q = zeros(ell,N);
for i = 1:ell
   zx = xdat.^pow(i,:);
   Q(i,:) = prod(zx,2);
end

% P matrix
pow = monpowers(2,maxDp);
m = size(pow,1); % number of nontrivial monomials in P
P = zeros(m,N);
for i = 1:m
   zy = ydat.^pow(i,:);
   P(i,:) = prod(zy,2)';
end

%% Create Koopman and Lie derivative matrix approximations

% Koopman approximation
K = P*pinv(Q);

% Symbolic variables
x = sdpvar(2,1); % 2D state variable
z = monolist(x,maxDp,0); % monomials that make up span(Dp)
w = monolist(x,maxDq,0); % monomials that make up span(Dq)
c = sdpvar(length(z),1); % coeffs to build Lyapunov function in span(Dp)

% Lie approximation
% ---> Lyapunov function = c.'*z in span(Dp)
L = c.'*(K*w) - c.'*z;

%% Identify Lyapunov function with SOS programming

% Inequality constraints posed relaxed to SOS conditions
cons = [sos(c.'*z - epsilon*dot(x,x)); sos(-L - epsilon*dot(x,x))];

% Set solver to MOSEK
opts = sdpsettings;
opts.solver = 'mosek';

% Objective function: minimize l^1 norm of coefficients
OBJ = sum(abs(c));

% Solve SOS problem
solvesos(cons,OBJ,opts,c)

% Print coefficients after removing those smaller than 10^-4 
c = clean(value(c), cleanVal)

%% Check identified Lyapunov function is a real Lyapunov function
% ---> This time we maximize epsilon and check if it is positive

% Initializations
yalmip clear
x = sdpvar(2,1); % 2D state variable
z = monolist(x,maxDp,0); % monomials that make up span(Dp)
sdpvar epsilon % epsilon is now a variable to be maximized

% Discovered Lyapunov function
v = c.'*z;

% True right-hand-side of map
f = [0.3*x(1); -x(1) + 0.5*x(2) + (7/18)*x(1)^2];

% True Lie derivative
Lie_exact = replace(v,x,f);

% Find Lyapunov function
solvesos([sos(v - epsilon*dot(x,x)); sos(- (Lie_exact - v) - epsilon*dot(x,x))],-epsilon)
fprintf('Maximized value of epsilon is: %f \n',value(epsilon))

if value(epsilon) > 0
   fprintf('We have discovered a true Lyapunov function from the data! \n') 
else
    fprintf('Something went wrong - we cannot verify that this is a Lyapunov function. \n')
end

















