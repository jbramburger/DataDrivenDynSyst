% -------------------------------------------------------------------------
% SINDy applied to Poincare and stroboscopic mapping discovery
%
% This script applies the SINDy method to discovering discrete-time
% iterative schemes that govern coarse-grained data from continuous-time 
% systems. Two examples are presented here:
%
% 1) A stroboscopic map of a periodically driven model for an RC circuit. 
% 2) A Poincare map of a planar model of the truncated Hopf normal form. 
%
% This script accompanies Section 3.2 of Data-Driven Methods for
% Dynamic Systems. 
%
% Author: Jason J. Bramburger
% -------------------------------------------------------------------------

% Clean workspace
clear all; close all; clc

%% Generate training data from RC circuit model

%Model parameters 
A = 1;
omega = 2*pi;

% Generate Trajectories RC circuit equation
dt = 0.01;
tspan = (0:20000-1)*dt;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,2));

%Generate More Trajectories
x0(1,:) = [0; 0]; %At least one solution
[~,xdat(1,:,:)]=ode45(@(t,x) RC(x,A,omega),tspan,x0(1,:),options);

numTraj = 3; %Number of trajectories
if numTraj >= 2
    for k = 2:numTraj
        
        x0(k,:) = [10*rand-5; 0]; %Initial conditions start in section 
        [~,xdat(k,:,:)]=ode45(@(t,x) RC(x,A,omega),tspan,x0(k,:),options);
        
    end
end

%% Aggregate stroboscopic section data

%Counting parameter
count = 1;

%Initialize
Psec(1) = xdat(1,1,1);
count = count + 1;

%Create Poincaré section data
for i = 1:numTraj
    for j = 1:length(xdat(i,:,1))-1 
        if (j == 1) && (i > 1) %Trajectories start in the section
            Psec(count) = xdat(i,j,1);
            count = count + 1; 
        elseif  (mod(xdat(i,j,2),2*pi/omega) >= 2*pi/omega-dt && mod(xdat(i,j+1,2),2*pi/omega) <= dt) 
            Psec(count) = xdat(i,j+1,1); %nth iterate
            PsecNext(count - 1) = xdat(i,j+1,1); %(n+1)st iterate
            count = count + 1;
        end
    end
    Psec = Psec(1:length(Psec)-1);
    count = count - 1;
end

% Create the recurrence data
xn = Psec;
xnp1 = PsecNext;

%% Create Theta matrix from monomials up to degree 5
 
Theta = ones([1,length(xn(1,:))]); % constant term
Theta(2,:) = xn; % linear term
Theta(3,:) = xn.^2; % quadratic term
Theta(4,:) = xn.^3; % cubic term
Theta(5,:) = xn.^4; % quartic term
Theta(6,:) = xn.^5; % quintic term

%% Find coefficients using pseudoinverse

Xi = xnp1*pinv(Theta);

%% Apply SINDy algorithm to obtain sparse model

% Sparsity parameter
lam = 0.01;

k = 1;
Xi_new = Xi;
while sum(sum(abs(Xi - Xi_new))) > 0  || k == 1 
    
    Xi = Xi_new;
    
    % find small coefficients
    smallinds = (abs(Xi) < lam); 
    
    % Threshold out small coefficients to eliminate from library
    Xi_new(smallinds) = 0;  
    
    % Identify the elements with large indices to remain in the library
    biginds = ~smallinds;
        
    % Find coefficients for reduced library
    Xi_new(biginds) = xnp1*pinv(Theta(biginds,:)); 
    
    k = k + 1;
end

%% Print model

% Clear command window
clc

% Monomials up to degree 5 in one variable
mons = ["" "x" "x^2" "x^3" "x^4" "x^5"];
XiString = string(abs(Xi_new));

fprintf('Discovered stroboscopic map using SINDy: \n')
fprintf('\n')

% Print xnp1 = f(xn) stroboscopic map model:
fprintf('f(x) = ')
bigcoeffs = abs(Xi_new(1,:)) > 1e-5; % chosen small just to weed out zero coefficients
for jnd = 1:length(bigcoeffs)
   if bigcoeffs(jnd) == 1
      
      % Print the model by excluding zeroed out terms 
      if Xi_new(1,jnd) < 0  
        fprintf('- %s ', strcat(XiString(jnd),mons(jnd)));
      else
        fprintf('+ %s ', strcat(XiString(jnd),mons(jnd)));
      end
      
   end
end
fprintf('\n')

%% Clean workspace for next experiment

clear all; close all; clc

%% Generate training data from planar Hopf normal form ODE

% model parameters 
om = 10*pi; % Hopf normal form period (make larger for better map result)

%Generate Trajectories Hopf normal form
m = 2; %Dimension of ODE
n = m-1; %Dimension of Poincare section
dt = 0.005;
tspan = (0:10000-1)*dt;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,m));

x0(1,:)=[0.0001; 0];  % Initial condition (Guarantee one starts near origin)
x0(2,:) = [0; 0]; %Trajectory at the fixed point

%Generate Trajectories
[~,xdat(1,:,:)]=ode45(@(t,x) Hopf(x,om),tspan,x0(1,:),options);
[~,xdat(2,:,:)]=ode45(@(t,x) Hopf(x,om),tspan,x0(2,:),options);

numTraj = 10;

if numTraj >= 3
    for k = 3:numTraj
        x0(k,:) = [5*rand; 0]; %Initial conditions start in section 
        [~,xdat(k,:,:)]=ode45(@(t,x) Hopf(x,om),tspan,x0(k,:),options);
    end
end

%% Aggregate Poincaré section data

%Counting parameter
count = 1;

%Initialize
Psec(1) = xdat(1,1,1);
count = count + 1;

%Create Poincare section data
for i = 1:numTraj
    for j = 1:length(xdat(1,:,1))-1 
        if  (xdat(i,j,2) < 0) && (xdat(i,j+1,2) >= 0) 
            Psec(count) = xdat(i,j+1,1); %nth iterate
            PsecNext(count - 1) = xdat(i,j+1,1); %(n+1)st iterate
            count = count + 1;
        end
    end
    Psec = Psec(1:length(Psec)-1);
    count = count - 1;
end

% Create the recurrence data
xn = Psec;
xnp1 = PsecNext;

%% Create Theta matrix from monomials up to degree 5 
%  --> One can increase or decrease the polynomial order with the parameter polyorder

Theta = ones([1,length(xn(1,:))]); % constant term

polyorder = 5; % change maximal polynomial order of monomials in library (print options below will not work if polyorder is not 5)

for p = 1:polyorder
    Theta(p+1,:) = xn.^p;
end

%% Find coefficients using pseudoinverse

Xi = xnp1*pinv(Theta);

%% Apply SINDy algorithm to obtain sparse model

% Sparsity parameter
lam = 0.01;

k = 1;
Xi_new = Xi;
while sum(sum(abs(Xi - Xi_new))) > 0  || k == 1 
    
    Xi = Xi_new;
    
    % find small coefficients
    smallinds = (abs(Xi) < lam); 
    
    % Threshold out small coefficients to eliminate from library
    Xi_new(smallinds) = 0;  
    
    % Identify the elements with large indices to remain in the library
    biginds = ~smallinds;
        
    % Find coefficients for reduced library
    Xi_new(biginds) = xnp1*pinv(Theta(biginds,:)); 
    
    k = k + 1;
end

%% Print model

% Clear command window
clc

% Monomials up to degree 5 in one variable
mons = ["" "x" "x^2" "x^3" "x^4" "x^5"];
XiString = string(abs(Xi_new));

fprintf('Discovered Poincare map using SINDy: \n')
fprintf('\n')

% Print xnp1 = f(xn) stroboscopic map model:
fprintf('f(x) = ')
bigcoeffs = abs(Xi_new(1,:)) > 1e-5; % chosen small just to weed out zero coefficients
for jnd = 1:length(bigcoeffs)
   if bigcoeffs(jnd) == 1
      
      % Print the model by excluding zeroed out terms 
      if Xi_new(1,jnd) < 0  
        fprintf('- %s ', strcat(XiString(jnd),mons(jnd)));
      else
        fprintf('+ %s ', strcat(XiString(jnd),mons(jnd)));
      end
      
   end
end
fprintf('\n')

%% RC circuit right-hand-side

function dx = RC(x,A,omega)

    dx = [A*sin(omega*x(2)) - x(1); 1];

end

%% Hopf Right-hand-side
function dx = Hopf(x,om)

    dx = [x(1) - om*x(2) - x(1)*(x(1)^2 + x(2)^2); om*x(1) + x(2) - x(2)*(x(1)^2 + x(2)^2)];

end

