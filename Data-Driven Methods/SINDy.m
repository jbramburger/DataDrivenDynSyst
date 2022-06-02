% -------------------------------------------------------------------------
% Sparse Identification of Nonlinear Dynamics (SINDy)
%
% This script applies the SINDy method of Proctor, Brunton, & Kutz (PNAS,
% 2016) to data collected from the Lorenz or Rossler chaotic systems. The
% script is broken into three parts: using a quadratic library, a cubic
% library, and then using the weak/integral SINDy formulation with a 
% quadratic library. 
%
% This script accompanies Section 2.2 of Computational Methods for
% Dynamical Systems. 
%
% To reproduce the noisy computations in the book load LorenzNoise.mat to
% get the noisy Lorenz data.
%
% Author: Jason J. Bramburger
% -------------------------------------------------------------------------

% Clean workspace
clear all; close all; clc

%% Generate synthetic training data from Lorenz or Rossler model

% Lorenz parameters 
sigma = 10;
beta = 8/3;
rho = 28;

% Rossler parameters
a = 0.1;
b = 0.1;
c = 18;

% Integration parameters
dt = 0.0001;
t = 0:dt:20;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,3));

% Simulate Lorenz system
x0 = [1;  2; rho - 1]; 
[t,xsol] = ode45(@(t,x) Lorenz(x,sigma,beta,rho),t,x0,options);

% Simulate Rossler system
% x0 = [0 -15 0];
% [t,xsol] = ode45(@(t,x) Rossler(x,a,b,c),t,x0,options);

% Add noise
var = 0.001; % noise variance
xsol = xsol + sqrt(var)*randn(size(xsol));
xsol = xsol'; % change dimensions to match theory

%% Load saved noisy data to reproduce result in the text

%load LorenzNoise.mat

%% Plot chaotic attractor

colormap(hot)
patch(xsol(1,:),xsol(2,:),xsol(3,:),t,'FaceColor','none','EdgeColor','interp','LineWidth',1) % rainbow coloured trajectory!
xlabel('$x$','Interpreter','Latex')
ylabel('$y$','Interpreter','Latex')
zlabel('$z$','Interpreter','Latex')
set(gca,'Fontsize',16)
axis tight
view(45,45)
grid on

%% Estimate derivatives

dxdt = (xsol(:,2:end) - xsol(:,1:end-1))/dt; % Y matrix
x = xsol(:,1:end-1); % X matrix

% Exact derivative for comparison
% for n = 1:length(t)-1
%    dxdt_exact(:,n) = Lorenz(x(;,n),sigma,beta,rho); 
% end

%% Create Theta matrix from monomials up to degree 2

Theta = ones([1,length(x(1,:))]); % constant term
Theta(2:4,:) = x; %linear term
Theta(5,:) = x(1,:).^2; % x^2 term
Theta(6,:) = x(1,:).*x(2,:); % xy term
Theta(7,:) = x(1,:).*x(3,:); % xz term
Theta(8,:) = x(2,:).^2; % y^2 term
Theta(9,:) = x(2,:).*x(3,:); % yz term
Theta(10,:) = x(3,:).^2; % z^2 term

%% Find coefficients using pseudoinverse

Xi = dxdt*pinv(Theta);

%% Apply SINDy algorithm to obtain sparse model

% Sparsity parameter
lam = 0.1;

k = 1;
Xi_new = Xi;
while sum(sum(abs(Xi - Xi_new))) > 0  || k == 1 
    
    Xi = Xi_new;
    
    % find small coefficients
    smallinds = (abs(Xi) < lam); 
    
    % Threshold out small coefficients to eliminate from library
    Xi_new(smallinds) = 0;  
    
    % Loop over all 3 variables of the Lorenz system
    for ind = 1:3 
        
        % Identify the elements with large indices to remain in the library
        biginds = ~smallinds(ind,:);
        
        % Find coefficients for reduced library
        Xi_new(ind,biginds) = dxdt(ind,:)*pinv(Theta(biginds,:)); 
    end
    
    k = k + 1;
end

%% Print model

% Clear command window
clc

% Monomials up to degree 2
mons2 = ["" "x" "y" "z" "x^2" "xy" "xz" "y^2" "yz" "z^2"];
XiString = string(abs(Xi_new));

fprintf('Discovered model using SINDy: \n')
fprintf('\n')

% Print dx/dt model:
fprintf('dx/dt = ')
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
fprintf('dy/dt = ')
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

% Print dy/dt model:
fprintf('dz/dt = ')
bigcoeffs = abs(Xi_new(3,:)) > 1e-5;
for jnd = 1:length(bigcoeffs)
   if bigcoeffs(jnd) == 1
      
      % Print the model by excluding zeroed out terms 
      if Xi_new(3,jnd) < 0  
        fprintf('- %s ', strcat(XiString(3,jnd),mons2(jnd)));
      else
        fprintf('+ %s ', strcat(XiString(3,jnd),mons2(jnd)));
      end
      
   end
end
fprintf('\n')

%% -------------------------------------------------------------------------
% Now add degree 3 monomials to the library too:
% -------------------------------------------------------------------------

%% Add cubic terms to library

Theta(11,:) = x(1,:).^3; % x^3 term
Theta(12,:) = x(1,:).*x(1,:).*x(2,:); % xxy term
Theta(13,:) = x(1,:).*x(2,:).*x(2,:); % xyy term
Theta(14,:) = x(1,:).*x(2,:).*x(3,:); % xyz term
Theta(15,:) = x(1,:).*x(3,:).*x(3,:); % xzz term
Theta(16,:) = x(2,:).*x(2,:).*x(2,:); % y^3 term
Theta(17,:) = x(2,:).*x(2,:).*x(3,:); % yyz term
Theta(18,:) = x(2,:).*x(3,:).*x(3,:); % yzz term
Theta(19,:) = x(3,:).*x(3,:).*x(3,:); % z^3 term

%% Find coefficients using pseudoinverse

Xi3 = dxdt*pinv(Theta);

%% Apply SINDy algorithm to obtain sparse model with cubic monomial terms 

% Sparsity parameter
lam = 0.1;

k = 1;
Xi3_new = Xi3;
while sum(sum(abs(Xi3 - Xi3_new))) > 0  || k == 1 
    
    Xi3 = Xi3_new;
    
    % find small coefficients
    smallinds = (abs(Xi3) < lam); 
    
    % Threshold out small coefficients to eliminate from library
    Xi3_new(smallinds) = 0;  
    
    % Loop over all 3 variables of the Lorenz system
    for ind = 1:3 
        
        % Identify the elements with large indices to remain in the library
        biginds = ~smallinds(ind,:);
        
        % Find coefficients for reduced library
        Xi3_new(ind,biginds) = dxdt(ind,:)*pinv(Theta(biginds,:)); 
    end
    
    k = k + 1;
end

%% Print model

% Clear command window
clc

% Monomials up to degree 3
mons3 = ["" "x" "y" "z" "x^2" "xy" "xz" "y^2" "yz" "z^2" "x^3" "x^2y" "xy^2" "xyz" "xz^2" "y^3" "y^2z" "yz^2" "z^3"];
Xi3String = string(abs(Xi3_new));

fprintf('Discovered model using SINDy: \n')
fprintf('\n')

% Print dx/dt model:
fprintf('dx/dt = ')
bigcoeffs = abs(Xi3_new(1,:)) > 1e-5;
for jnd = 1:length(bigcoeffs)
   if bigcoeffs(jnd) == 1
      
      % Print the model by excluding zeroed out terms 
      if Xi_new(1,jnd) < 0  
        fprintf('- %s ', strcat(Xi3String(1,jnd),mons3(jnd)));
      else
        fprintf('+ %s', strcat(Xi3String(1,jnd),mons3(jnd)));
      end
      
   end
end
fprintf('\n')

% Print dy/dt model:
fprintf('dy/dt = ')
bigcoeffs = abs(Xi3_new(2,:)) > 1e-5;
for jnd = 1:length(bigcoeffs)
   if bigcoeffs(jnd) == 1
      
      % Print the model by excluding zeroed out terms 
      if Xi3_new(2,jnd) < 0  
        fprintf('- %s ', strcat(Xi3String(2,jnd),mons3(jnd)));
      else
        fprintf('+ %s ', strcat(Xi3String(2,jnd),mons3(jnd)));
      end
      
   end
end
fprintf('\n')

% Print dy/dt model:
fprintf('dz/dt = ')
bigcoeffs = abs(Xi3_new(3,:)) > 1e-5;
for jnd = 1:length(bigcoeffs)
   if bigcoeffs(jnd) == 1
      
      % Print the model by excluding zeroed out terms 
      if Xi3_new(3,jnd) < 0  
        fprintf('- %s ', strcat(Xi3String(3,jnd),mons3(jnd)));
      else
        fprintf('+ %s ', strcat(Xi3String(3,jnd),mons3(jnd)));
      end
      
   end
end
fprintf('\n')

%% -------------------------------------------------------------------------
% Now use integrals instead of derivatives to counter noise:
% -------------------------------------------------------------------------

%% Set up matrices

x = xsol(:,1:end);
Y = xsol(:,2:end) - xsol(:,1).*ones(size(xsol(:,2:end)));
ThetaInt = zeros(size(Theta(1:10,:))); % Just keep quadratic terms

% Initialize ThetaInt
ThetaInt(1:10,1) = dt*Theta(1:10,1);

for  n = 2:length(Theta(1,:))
    ThetaInt(1:10,n) = ThetaInt(1:10,n-1) + dt*Theta(1:10,n);
end

%% Initialize SINDy with least-squares solution

XiInt = Y*pinv(ThetaInt);

%% Apply SINDy algorithm to obtain sparse model with quadratic monomial integral terms 

% Sparsity parameter
lam = 0.1;

k = 1;
XiInt_new = XiInt;
while sum(sum(abs(XiInt - XiInt_new))) > 0  || k == 1 
    
    XiInt = XiInt_new;
    
    % find small coefficients
    smallinds = (abs(XiInt) < lam); 
    
    % Threshold out small coefficients to eliminate from library
    XiInt_new(smallinds) = 0;  
    
    % Loop over all 3 variables of the Lorenz system
    for ind = 1:3 
        
        % Identify the elements with large indices to remain in the library
        biginds = ~smallinds(ind,:);
        
        % Find coefficients for reduced library
        XiInt_new(ind,biginds) = Y(ind,:)*pinv(ThetaInt(biginds,:)); 
    end
    
    k = k + 1;
end 

%% Print Integral SINDy model

% Clear command window
clc

% Monomials up to degree 2
XiIntString = string(abs(XiInt_new));

fprintf('Discovered model using SINDy: \n')
fprintf('\n')

% Print dx/dt model:
fprintf('dx/dt = ')
bigcoeffs = abs(XiInt_new(1,:)) > 1e-5;
for jnd = 1:length(bigcoeffs)
   if bigcoeffs(jnd) == 1
      
      % Print the model by excluding zeroed out terms 
      if XiInt_new(1,jnd) < 0  
        fprintf('- %s ', strcat(XiIntString(1,jnd),mons2(jnd)));
      else
        fprintf('+ %s', strcat(XiIntString(1,jnd),mons2(jnd)));
      end
      
   end
end
fprintf('\n')

% Print dy/dt model:
fprintf('dy/dt = ')
bigcoeffs = abs(XiInt_new(2,:)) > 1e-5;
for jnd = 1:length(bigcoeffs)
   if bigcoeffs(jnd) == 1
      
      % Print the model by excluding zeroed out terms 
      if XiInt_new(2,jnd) < 0  
        fprintf('- %s ', strcat(XiIntString(2,jnd),mons2(jnd)));
      else
        fprintf('+ %s ', strcat(XiIntString(2,jnd),mons2(jnd)));
      end
      
   end
end
fprintf('\n')

% Print dy/dt model:
fprintf('dz/dt = ')
bigcoeffs = abs(XiInt_new(3,:)) > 1e-5;
for jnd = 1:length(bigcoeffs)
   if bigcoeffs(jnd) == 1
      
      % Print the model by excluding zeroed out terms 
      if XiInt_new(3,jnd) < 0  
        fprintf('- %s ', strcat(XiIntString(3,jnd),mons2(jnd)));
      else
        fprintf('+ %s ', strcat(XiIntString(3,jnd),mons2(jnd)));
      end
      
   end
end
fprintf('\n')

%% Lorenz right-hand-side

function dx = Lorenz(x,sigma,beta,rho)


    dx = [sigma*(x(2) - x(1)); x(1)*(rho - x(3)) - x(2); x(1)*x(2) - beta*x(3)];
    
end

%% Rossler right-hand-side

function dx = Rossler(x,a,b,c)

    dx = [-x(2) - x(3); x(1) + a*x(2); b + x(3)*(x(1) - c)];

end