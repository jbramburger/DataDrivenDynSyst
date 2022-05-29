% -------------------------------------------------------------------------
% Using SINDy to perform data-driven averaging
%
% This script applies the method of data-driven averaging to multiscale
% time series data. The method relies on the assumption that the fast scale
% is periodic and leverages the SINDy method to discover the slow scale
% drift in the data. Two examples are provided in this script. The first is
% the periodically forced logistic model for which results can be compared
% with analytical results from the theory of averaging for dynamical
% systems. The second applies the method to the motions of Jupiter or
% Saturn in their orbital plane about the sun.
%
% This script accompanies Section 3.3 of Data-Driven Methods for Dynamic 
% Systems. 
%
% Jupiter and Saturn orbits are read in from the jupiter_data.mat and
% saturn_data.mat files, respectively.
%
% Author: Jason J. Bramburger
% -------------------------------------------------------------------------

% Clean workspace
clear all; close all; clc

%% Simulate the periodically forced logistic model

%Model parameters
eps = 0.001; % epsilon value
T = 2*pi; % period of forcing terms

%Inializations for generating trajectories
m = 2; %Dimension of ODE (including time)
dt = 0.005;
tspan = (0:500000-1)*dt; % need long timescale to see slow dynamics
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,m));

% Generate trajectory
x0 = [4*rand; 0];  
[~,xdat]=ode45(@(t,x) Logistic(x,eps,T),tspan,x0,options);

%% Plot results
plot(xdat(:,2),xdat(:,1),'b','LineWidth',1)
xlabel('t')
ylabel('x(t)')
set(gca,'Fontsize',16)

%% Collect coarse-grained data

% Counting parameter
count = 1;

% Known equilibrium entry
xn(1) = 0;
xnp1(1) = 0;

% Collecting mapping data
for j = 1:length(xdat(:,1))-1 
    if  (mod(xdat(j,2),1) >= 1 - 2*dt && mod(xdat(j+1,2),1) <= 2*dt) %&& j >= length(xdat(i,:,1))/2) 
        temp(count) = xdat(j+1,1); %nth iterate
        count = count + 1;
    end
end
xn = [xn; temp(1:length(temp)-1)']';
xnp1 = [xnp1; temp(2:length(temp))']';

%% Apply SINDy method for map discovery

Theta = ones([1,length(xn)]); % constant term
Theta(2,:) = xn; % linear term
Theta(3,:) = xn.^2; % quadratic term
Theta(4,:) = xn.^3; % cubic term
Theta(5,:) = xn.^4; % quartic term
Theta(6,:) = xn.^5; % quintic term

% Find coefficients using pseudoinverse
Xi = xnp1*pinv(Theta);

% Sparsity parameter
lam = eps^2;

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

% Clear command window
clc

% Monomials up to degree 5 in one variable
mons = ["" "x" "x^2" "x^3" "x^4" "x^5"];
XiString = string(abs(Xi_new));

fprintf('Discovered discrete-time flow map using SINDy: \n')
fprintf('\n')

% Print xnp1 = f(xn) flow map model:
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

%% Determine slow ODE

Xi_slow = Xi_new;
Xi_slow(2) = Xi_slow(2) - 1; % subtract off identity
Xi_slow = Xi_slow/eps;
XiString = string(abs(Xi_slow));

fprintf('Discovered continuous-time slow-scale ODE from flow-map: \n')
fprintf('\n')

% Print dx/dt = g(x) slow drift model:
fprintf('g(x) = ')
bigcoeffs = abs(Xi_slow(1,:)) > 1e-5; % chosen small just to weed out zero coefficients
for jnd = 1:length(bigcoeffs)
   if bigcoeffs(jnd) == 1
      
      % Print the model by excluding zeroed out terms 
      if Xi_slow(1,jnd) < 0  
        fprintf('- %s ', strcat(XiString(jnd),mons(jnd)));
      else
        fprintf('+ %s ', strcat(XiString(jnd),mons(jnd)));
      end
      
   end
end
fprintf('\n')

%% Clean workspace for next experiment

clear all; close all; clc

%% Load Jupiter/Saturn data for eccentricity discovery

format long

% Load Jupiter data
load jupiter_data;
data = jupiter_data;
T = 11.85; % Jupiter fast period

% Load Saturn data <----- Uncomment below for Saturn data
% load saturn_data;
% data = saturn_data;
% T = 29.5; % Saturn fast period

% Initializations
n = 2; % number of components to F and G functions
l = length(tspan);
dt = tspan(2) - tspan(1);
sample = [1 2]; % restrict to orbital plane
count = 1;

%% Coarsened Data

for j = 1:l-1
   if mod(tspan(j+1),T) == 0 
        xn(1:n,count) = data(sample,j);
        count = count+1;
   end
end

cut = 1500; % Truncate size of input data 
xnp1 = xn(:,2:cut);
xn = xn(:,1:cut-1);

%% Apply SINDy Method to discover flow-map (F)

% Generate cubic library
Theta = ones([1,length(xn(1,:))]); % constant term
Theta(2:3,:) = xn; % linear terms
Theta(4,:) = xn(1,:).^2; % x1^2
Theta(5,:) = xn(1,:).*xn(2,:); % x1*x2
Theta(6,:) = xn(2,:).^2; % x2^2
Theta(7,:) = xn(1,:).^3; % x1^3
Theta(8,:) = xn(2,:).*xn(1,:).^2; % x2*x1^2
Theta(9,:) = xn(1,:).*xn(2,:).^2; % x1*x2^2
Theta(10,:) = xn(2,:).^3; % x2^3

% Find coefficients using pseudoinverse
Xi = xnp1*pinv(Theta);

% Sparsity parameter
lam = 1e-3;

k = 1;
Xi_new = Xi;
while sum(sum(abs(Xi - Xi_new))) > 0  || k == 1 
    
    Xi = Xi_new;
    
    % find small coefficients
    smallinds = (abs(Xi) < lam); 
    
    % Threshold out small coefficients to eliminate from library
    Xi_new(smallinds) = 0;  
    
    % Loop over both x and y variables in the orbital plane
    for ind = 1:n 
        
        % Identify the elements with large indices to remain in the library
        biginds = ~smallinds(ind,:);
        
        % Find coefficients for reduced library
        Xi_new(ind,biginds) = xnp1(ind,:)*pinv(Theta(biginds,:)); 
    end
    
    k = k + 1;
end

% Clear command window
clc

% Monomials up to degree 2
mons2 = ["" "x" "y" "x^2" "xy" "y^2" "x^3" "x^2y" "xy^2" "y^3"];
XiString = string(abs(Xi_new));

fprintf('Discovered model using SINDy: \n')
fprintf('\n')

% Print dx/dt model:
fprintf('F1(x,y) = ')
bigcoeffs = abs(Xi_new(1,:)) > 1e-5; % chosen small just to weed out zero coefficients
for jnd = 1:length(bigcoeffs)
   if bigcoeffs(jnd) == 1
      
      % Print the model by excluding zeroed out terms 
      if Xi_new(1,jnd) < 0  
        fprintf('- %s ', strcat(XiString(1,jnd),mons2(jnd)));
      else
        fprintf('+ %s', strcat(XiString(1,jnd),mons2(jnd)));
      end
      
   end
end
fprintf('\n')

% Print dy/dt model:
fprintf('F2(x,y) = ')
bigcoeffs = abs(Xi_new(2,:)) > 1e-5;
for jnd = 1:length(bigcoeffs)
   if bigcoeffs(jnd) == 1
      
      % Print the model by excluding zeroed out terms 
      if Xi_new(2,jnd) < 0  
        fprintf('- %s ', strcat(XiString(2,jnd),mons2(jnd)));
      else
        fprintf('+ %s ', strcat(XiString(2,jnd),mons2(jnd)));
      end
      
   end
end
fprintf('\n')

%% Determine eccentricity ODE (G)

eps = T/46800; % scale separation parameter
Xi_slow = Xi_new;
Xi_slow(1,2) = 0; % subtract off identity (rounded)
Xi_slow(2,3) = 0; % subtract off identity
Xi_slow = Xi_slow/(T*eps);
XiString = string(abs(Xi_slow));

fprintf('Discovered continuous-time slow-scale ODE from flow-map: \n')
fprintf('\n')

% Print slow drift model:
% First component:
fprintf('G1(x,y) = ')
bigcoeffs = abs(Xi_slow(1,:)) > 1e-5; % chosen small just to weed out zero coefficients
for jnd = 1:length(bigcoeffs)
   if bigcoeffs(jnd) == 1
      
      % Print the model by excluding zeroed out terms 
      if Xi_slow(1,jnd) < 0  
        fprintf('- %s ', strcat(XiString(1,jnd),mons2(jnd)));
      else
        fprintf('+ %s', strcat(XiString(1,jnd),mons2(jnd)));
      end
      
   end
end
fprintf('\n')

% Second component:
fprintf('G2(x,y) = ')
bigcoeffs = abs(Xi_slow(2,:)) > 1e-5;
for jnd = 1:length(bigcoeffs)
   if bigcoeffs(jnd) == 1
      
      % Print the model by excluding zeroed out terms 
      if Xi_slow(2,jnd) < 0  
        fprintf('- %s ', strcat(XiString(2,jnd),mons2(jnd)));
      else
        fprintf('+ %s ', strcat(XiString(2,jnd),mons2(jnd)));
      end
      
   end
end
fprintf('\n')

%% Logistic right-hand-side

function dx = Logistic(x,eps,T)

    dx = [eps*(x(1)*(1- x(1) + sin(T*x(2)))); 1];

end






