% -------------------------------------------------------------------------
% Identifying Lyapunov Functions with Semidefinite Programming
%
% The goal of this script is to use semidefinite programming to identify a
% Lyapunov function for the planar Moore-Greitzer system
%      x1' = -x2 - 1.5x1^2 - 0.5x1^3, 
%      x2' = 3x1 - x2.
% The script expands the unknown Lyapunov function as a polynomial up to a
% given degree and tunes the coefficients in order to satisfy both of the
% constraints:
%      V - x1^2 - x2^2 is SOS
%      -LV - x1^2 - x1^2 is SOS
% where LV is the Lie derivative associated to the Moore-Greitzer system.
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

%% Discovering polynomial Lyapunov function

% Degree of Lyapunov function
degV = 6;

% Symbolic variables
x = sdpvar(2,1);

% Moore-Greitzer right-hand-side
f = [-x(2) - 1.5*x(1)^2 - 0.5*x(1)^3; 3*x(1) - x(2)];

% Lyapunov function, gradient, and Lie derivative
[v,vc,vm] = polynomial(x,degV,1); 
Dv = jacobian(v,x);        
Lie = dot(Dv,f);

% Find Lyapunov function
cons = [sos(v - dot(x,x)); sos(-Lie - dot(x,x))]; % SOS constraints to be satisfied 
obj = sum(abs(vc)); % objective function: minimize l^1 norm of coefficients 
sol = solvesos(cons,obj,[],vc);

%% Print Results

% Clean window
clc

% Print if solution was solved
if sol.problem == 0 % <-- solution is solved successfully
    
    % Remove out small values from coefficients
    coeffs = clean(value(vc),1e-5);
    
    % Display Lyapunov function
    fprintf('Lyapunov function identified! \n')
    fprintf('\n')
    fprintf('V(x(1),x(2)) = ') 
    sdisplay(coeffs.'*vm)
else
    % If no Lyapunov function, tell the user
    fprintf('Sorry, unable to identify a Lyapunov function \n')
end

%% Plot Lyapunov function

% Spatial coordinates
N = 1000;
x2 = linspace(-5,5,N);
y2 = linspace(-5,5,N);
[X,Y] = meshgrid(x2,y2);

% Build Lyapunov function
v_num = dot(value(vc),vm);
[pow, coeff] = getexponentbase(v_num,x);
lyap = zeros(size(X));
for i = 1:numel(lyap)
    lyap(i) = lyap(i) + coeff*( X(i).^pow(:,1) .* Y(i).^pow(:,2) );
end

% Surface plot
surf(X,Y,lyap)
shading interp
xlabel('$x_1$','Interpreter','Latex')
ylabel('$x_2$','Interpreter','Latex')
zlabel('$V(x_1,x_2)$','Interpreter','Latex')
set(gca,'FontSize',16,'Xlim',[-5 5],'Ylim',[-5 5])


