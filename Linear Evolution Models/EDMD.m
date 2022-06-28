% -------------------------------------------------------------------------
% Extended Dynamic Mode Decomposition (EDMD)
%
% This code applies the Extended DMD (EDMD) approach to a simple nonlinear
% discrete dynamical system for which the Koopman eigenfunctions are known.
% Below we compare results using a dictionary of monomials up to degree 2
% and up to degree 3. Training data is obtained by drawing N (number of
% snapshots) pairs of points from the uniform distribution on [-1,1]^2 and
% tracking their image under the dynamical system. The system depends on
% two parameters, lam1 and lam2, which can be varied to produce similar
% results.
%
% This script accompanies Section 2.3 of Computational Methods for
% Dynamical Systems. 
%
% Author: Jason J. Bramburger
% -------------------------------------------------------------------------

% Clean workspace
clear all; close all; clc

%% Generate snapshot data

N = 10000; % numer of snapshots
X = 2*rand(2,N) - 1;
Y = zeros(2,N);

% System parameters
lam1 = 1.5;
lam2 = 0.5;

% Populate Y matrix
for n = 1:N
    Y(1,n) = lam1*X(1,n);
    Y(2,n) = lam2*(X(2,n) - X(1,n)^2);
end

%% Create matrices of observables Psi(X) and Psi(Y)
% Start with observable functions of monomials up to degree 2

% First two rows are degree 1 variables (x_1 and x_2)
PsiX = X;
PsiY = Y;

% x_1^2 observable 
PsiX(3,:) = X(1,:).^2; 
PsiY(3,:) = Y(1,:).^2; 

% x_1*x_2 observable 
PsiX(4,:) = X(1,:).*X(2,:); 
PsiY(4,:) = Y(1,:).*Y(2,:); 

% x_2^2 observable 
PsiX(5,:) = X(2,:).^2; 
PsiY(5,:) = Y(2,:).^2;

%% DMD matrix

A = PsiY*pinv(PsiX);
[V, D, W] = eig(A); % compute eigenvalues + left (W) & right (V) eigenvectors
mu = diag(D); % extract eigenvalues

%% Koopman Eigenfunctions

% Clean window
clc

% Print approximate eigenfunctions
for p = 1:5
   fprintf('Eigenvalue: %4.4f \n',mu(p))
   fprintf('Associated (approximate) Koopman eigenfunction:\n')
   fprintf('%4.4f*x1 + %4.4f*x2 + %4.4f*x1^2 + %4.4f*x1*x2 + %4.4f*x2^2\n', W(:,p))
   fprintf('\n')
end

%% Add in cubic observables

% x_1^3 observable 
PsiX(6,:) = X(1,:).^3; 
PsiY(6,:) = Y(1,:).^3; 

% x_1^2*x_2 observable 
PsiX(7,:) = X(1,:).*X(1,:).*X(2,:); 
PsiY(7,:) = Y(1,:).*Y(1,:).*Y(2,:); 

% x_1*x_2^2 observable 
PsiX(8,:) = X(1,:).*X(2,:).*X(2,:); 
PsiY(8,:) = Y(1,:).*Y(2,:).*Y(2,:); 

% x_2^3 observable 
PsiX(9,:) = X(2,:).^3; 
PsiY(9,:) = Y(2,:).^3; 

%% DMD matrix

A3 = PsiY*pinv(PsiX);
[V3, D3, W3] = eig(A3); % compute eigenvalues + left (W) & right (V) eigenvectors
mu3 = diag(D3); % extract eigenvalues

%% Koopman Eigenfunctions

% Clean window
clc

% Print approximate eigenfunctions
for p = 1:9
   fprintf('Eigenvalue: %4.5f \n',mu3(p))
   fprintf('Associated (approximate) Koopman eigenfunction:\n')
   fprintf('%4.5f*x1 + %4.5f*x2 + %4.5f*x1^2 + %4.5f*x1*x2 + %4.5f*x2^2 + %4.5f*x1^3 + %4.5f*x1^2*x2 + %4.5f*x1*x2^2 + %4.5f*x2^3\n', W3(:,p))
   fprintf('\n')
end







