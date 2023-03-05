% -------------------------------------------------------------------------
% Time Average Bounds on the Stochastic Logistic Map
%
% The goal of the script is to use data to approximate the Lie derivative 
% associated to the iterations of the stochastic logistic map to bound 
% long-time (expected) average the state-space observable, Phi(x) = x.
%
% Instead of the usual monomial basis, we use a Chebyshev basis to improve 
% numerical conditioning and accuracy. This requires one to use the ChebFun 
% package, which can be downloaded at: https://www.chebfun.org/download/
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
clear; 
close all; 
clc
yalmip clear
format long

%% Method Parameters 
% al = max degree of Dp dictionary of obserables
% degF = assumed degree of underlying map (2 in the case of logistic)
% m = max degree of Dq dictionary of obserables
al = 2;
degF = 2; 
beta = degF*al;
N = 1e6;
TOL = 1/sqrt(N);    % thresholding tolerance for the EDMD matrix

%% Generate synthetic data

% Chebfun objects -- note the data is in [0,1]
T1 = chebpoly(0:al,[0,1]);
T2 = chebpoly(0:beta,[0,1]);

% Generate data
x = zeros(N,1);
x(1) = rand;
for n = 2:N
    % parameter is drawn from uniform distribution on [0,4] at each
    % iteration
    x(n) = 4*rand*x(n-1)*(1 - x(n-1)); 
end

% Load data to replicate results from manuscript
% load('stoch_logistic_traj_N1e+04.mat')
% N = length(x);

% Compute emperical average 
avg = sum(x)/N;

% Koopman matrix
Phi = T2( x(1:N-1) )'; 
Psi = T1( x(2:N) )'; 
K = edmd_with_thresholding(Phi,Psi,TOL);

%% Compute upper bounds from data

% Matrices and sdp variables to SOS problem
yalmip clear
Q = sdpvar(beta/2+1,beta/2+1);
P = sdpvar(beta/2,beta/2);
c = sdpvar(al+1,1);
sdpvar B
p = K'*c;                   % coeffs of Kv
p(1:al+1) = p(1:al+1) - c;    % coeffs of Lv = Kv - v
p(1:2) = p(1:2) + 0.5;            % coeffs of Lv + x
p(1) = p(1) - B;            % coeffs of Lv + x - B
p = -p;                     % coeffs of the expression to be SOS
At = chebsdp_1d(beta/2);
Bt = chebsdp_1d_locball(beta/2);

% SOS problem setup and solving
CNSTR = [Q>=0, P>=0, At.'*Q(:)+Bt.'*P(:)==p];    % SDP constraints
opts = sdpsettings('dualize',1); % problem in primal standard form
optimize(CNSTR,B,opts)
Bu = value(B);

%% Compute lower bounds from data

% Matrices and sdp variables to SOS problem
yalmip clear
Q = sdpvar(beta/2+1,beta/2+1);
P = sdpvar(beta/2,beta/2);
c = sdpvar(al+1,1);
sdpvar B
p = K'*c;                   % coeffs of Kv
p(1:al+1) = p(1:al+1) - c;    % coeffs of Lv = Kv - v
p(1:2) = p(1:2) + 0.5;            % coeffs of Lv + x
p(1) = p(1) - B;            % coeffs of Lv + x - B
At = chebsdp_1d(beta/2);
Bt = chebsdp_1d_locball(beta/2);

% SOS problem setup and solving
CNSTR = [Q>=0, P>=0, At.'*Q(:)+Bt.'*P(:)==p];    % SDP constraints
opts = sdpsettings('dualize',1); % problem in primal standard form
optimize(CNSTR,-B,opts)
Bl = value(B);

%% Display results

% clean window
clc

% Print bounds
fprintf('Numerically computed upper bound: %f \n',value(Bu))
fprintf('Emperical average computed from data: %f \n',avg)
fprintf('Numerically computed lower bound: %f \n',value(Bl))
fprintf('\n')

% Symbolic koopman from data
TOL = 1e-10;
fprintf('EDMD Koopman:\n')
for j = 0:al
    line = sprintf('K( T_%i(x) ) = ', j);
    if abs(K(j+1,1))>1e-10
       next = sprintf(' %+8.6f', K(j+1,1));
       line = [line, next];
    end
    for i = 1:beta
        if K(j+1,i+1) >= TOL
            next = sprintf(' + %8.6f T_%i(x)',K(j+1,i+1),i);
        elseif K(j+1,i+1) <= -TOL
            next = sprintf(' - %8.6f T_%i(x)',-K(j+1,i+1),i);
        else
            next = [];
        end
        line = [line, next];
    end
    fprintf([line, '\n'])
end
fprintf('\n')

% Exact koopman
lam = chebfun2(@(lam,x) lam,[0,4,0,1]);
x = chebfun2(@(lam,x) x,[0,4,0,1]);
fprintf('Exact Koopman:\n')
for j = 0:al
    T = chebpoly(j, [0,1]);
    KT = sum(T(lam.*x.*(1-x)), 2); % Integrate chebfun over l
    c = chebcoeffs(KT)./4;
    line = sprintf('K( T_%i(x) ) = ', j);
    if abs(c(1))>1e-10
       next = sprintf(' %+8.6f', c(1));
       line = [line, next];
    end
    for i = 2:length(c)
        if c(i) >= TOL
            next = sprintf(' + %8.6f T_%i(x)',c(i),i-1);
        elseif c(i) <= -TOL
            next = sprintf(' - %8.6f T_%i(x)',-c(i),i-1);
        else
            next = [];
        end
        line = [line, next];
    end
    fprintf([line, '\n'])
end




