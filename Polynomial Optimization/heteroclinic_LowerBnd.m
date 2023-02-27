% -------------------------------------------------------------------------
% Lower Bounds on the Existence of a Heteroclinic Orbit
%
% The goal of this script is to use semidefinite programming to provide
% lower bounds on the range of mu values for which a heteroclinic 
% connection exists in the system
%      x1' = x2, 
%      x2' = -mu*x2 - x1*(1-x1).
% The script searches for a desired auxiliary function using SOS 
% relaxations and localizations in phase space.   
%
% To run this notebook one requires paths to the freely available software
% packages YALMIP and MOSEK. YALMIP can be downloaded at:
%           https://yalmip.github.io/download/
% and MOSEK can be downloaded at:
%           https://www.mosek.com/downloads/
%
% This script accompanies Section 4.2 of Data-Driven Methods for
% Dynamic Systems. 
%
% Author: Jason J. Bramburger
% -------------------------------------------------------------------------

% Clean workspace
clear all; 
close all;
yalmip clear
clc

%% Bisection Method for Finding Upper Bound

format long

% ODE mu parameter
m = 2; 

% Bounding method parameters
degV = 20; % Degree of auxiliary function
lambda = 10^3;

%Bisection Method
muleft = 0;
muright = 2;
while abs(muright - muleft) >= 1e-5

    mumid = 0.5*(muleft + muright);
    flag = volume(mumid,m,degV,lambda);
    
    if flag == 0
       muleft = mumid; 
    else
       muright = mumid;
    end
    
    mumid
end

%Print Results
fprintf('The lower bound on the minimum speed for m = %d is %f found using degree %d polynomials.\n',m,mumid,degV)

%%
function flag = volume(mu,m,d,lambda)

% SDP variables
sdpvar u v

% Epsilon value
eps = 1e-4;

% Auxiliary function
[V, cV] = polynomial([u v], d);

% S Procedure
d2 = d+m;
[s1, c1] = polynomial([u v], d2);
[s2, c2] = polynomial([u v], d2);
[s3, c3] = polynomial(u, d2);

% Derivatives
dVdu = jacobian(V,u);
dVdv = jacobian(V,v);

% Replacements
Vv0 = replace(V,v,0);

%Constraints
cons = [];
cons = [cons, replace(V,[u v], [0 -eps]) == 0, replace(V, [u v], [1 0]) == 0];
cons = [cons, sos(lambda*(dVdu*v - mu*dVdv*v - dVdv*(1-u)*u^m) - V - u*(1-u)*s1 + v*s2)]; 
cons = [cons, sos(-Vv0 - eps*(1-u) - u*(1-u)*s3)];
cons = [cons, sos(s1), sos(s2), sos(s3)];

%SOS Solver
ops = sdpsettings('solver','mosek','verbose',0,'cachesolvers',1);
sol = solvesos(cons,[],ops,[cV;c1; c2; c3]);

%Return whether solvesos failed or succeeded
flag = sol.problem;

end














