% -------------------------------------------------------------------------
% Learning Parameter-Dependent Poincare Maps in the Sprott System
%
% We use the third-order chaotic ordinary differential equation:
%
%        x''' + mu*x'' - (x')^2 + x = 0,
%
% where mu is the bifurcation parameter. The user is directed to the paper
% "Simplest dissipative chaotic flow" by J.C. Sprott (Phys. Lett. A, 228,
% pp. 271-274, 1997) for a complete discussion of the system dynamics and 
% its bifurcations. 
%
% TBD...
%
% This script accompanies Section 3.3 of Data-Driven Methods for
% Dynamic Systems. 
%
% For applications of this method to other systems see the repository:
%      https://github.com/jbramburger/Stabilizing_UPOs
%
% Author: Jason J. Bramburger
% -------------------------------------------------------------------------

% Clean workspace
clear all; close all; clc

%% Add path to YALMIP (assumed to be in the folder 'YALMIP-master')

addpath(genpath('YALMIP-master'))

%% Gather data from the Sprott system

% Bifurcation parameter 
mustar = 2.06;
%mu = 2.08:0.001:2.12; % focal parameter = mustar = 2.10 (period 2 attractor)
mu = 2.04:0.001:2.08; % focal parameter = mustar = 2.06 (period 8 attractor)

%ODE generation parameters
m = 3; % Dimension of ODE 
n = m; % Dimension of Poincare section (2D + 1D parameter)
dt = 0.01;
tspan = (0:1000000-1)*dt;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,m));

%Generate Trajectories
for k = 1:length(mu)
   
    x0 = [6; 0; -2.5]; 
   [~,xdat(k,:,:)]=ode45(@(t,x) Sprott(x,mu(k)),tspan,x0,options); 
   
end

%% Poincare section data (y = 0 with y' < 0)

%Counting parameter
count = 1;

Psec = [];
PsecNext = [];

%Create Poincare section data
for i = 1:length(mu)
   temp = [];
    for j = 2:length(xdat(i,:,1))-1 
        if  (xdat(i,j,2) >= 0 && xdat(i,j+1,2) <= 0) % section condition 
            temp(count,:) = [xdat(i,j,1) xdat(i,j,3) mu(i)]; %nth iterate
            count = count + 1;
        end
    end
    Psec = [Psec; temp(1:length(temp)-1,:)];
    PsecNext = [PsecNext; temp(2:length(temp),:)];
   	count = 1; 
end

% Shift data
shift = mean(Psec);
Psec = Psec - shift;
PsecNext = PsecNext - shift;

%% Create Theta matrix

polyorder = 3;
pow = monpowers(3,polyorder);
M = size(pow,1); % number of nontrivial monomials
N = size(Psec,1);
Theta = zeros(N,M);
for i = 1:M
    z = Psec.^pow(i,:);
    Theta(:,i) = prod(z,2);
end

%% SINDy procedure

Xi = pinv(Theta)*PsecNext;

% Sparsity parameter
lam = 0.01;

% Sparsity promotion using SINDy
k = 1;
Xi_new = Xi;
while sum(sum(abs(Xi - Xi_new))) > 0  || k == 1 
    
    Xi = Xi_new;
    
    % find small coefficients
    smallinds = (abs(Xi) < lam); 
    
    % Threshold out small coefficients to eliminate from library
    Xi_new(smallinds) = 0;  
    
    % Loop over all 3 variables in the Poincare section
    for ind = 1:n 
        
        % Identify the elements with large indices to remain in the library
        biginds = ~smallinds(:,ind);
        
        % Find coefficients for reduced library
        Xi_new(biginds,ind) = pinv(Theta(:,biginds))*PsecNext(:,ind); 
    end
    
    k = k + 1;
end

%% Build SINDy map with symbolic variables

syms x z p

% xnp1 = f(xn,zn) and znp1 = g(xn,zn)
f(x,z,p) = dot(prod([x z p].^pow,2),Xi_new(:,1));
g(x,z,p) = dot(prod([x z p].^pow,2),Xi_new(:,2));

%% Identify fixed point 

% mustar value in shifted coordinates of the SINDy mapping
mu_s = mustar - shift(3);

% Replace parameter p with mu_s
fs(x,z) = f(x,z,mu_s);
gs(x,z) = g(x,z,mu_s);

% Identify (shifted) fixed point
fixed = vpasolve([fs(x,z) == x, gs(x,z) == z], [x, z]);
x1 = double(fixed.x);
z1 = double(fixed.z);

% Symbolic derivatives
dfx = diff(f,x);
dfz = diff(f,z);
dfp = diff(f,p);
dgx = diff(g,x);
dgz = diff(g,z);
dgp = diff(g,p);

% Jacobian matrices
A = double([dfx(x1,z1,mu_s), dfz(x1,z1,mu_s); dgx(x1,z1,mu_s), dgz(x1,z1,mu_s)]);
B = double([dfp(x1,z1,mu_s); dgp(x1,z1,mu_s)]);

%% Identify control matrix

% Epsilon value for numerical strict inequalities
eps = 1e-3;
I = eps*eye(4);

% SDP variables
Q = sdpvar(2,2);
Y = sdpvar(1,2,'full');

% Semidefinite constraints and optimization
M = [Q Q*A' + Y'*B'; A*Q+B*Y Q];
lmi = [M - I >= 0];
optimize(lmi)

% Control matrix for fixed points
C1 = value(Y)/(value(Q)); 

%% Controlled trajectory (fixed point)

% Shift variables back
x1 = double(fixed.x) + shift(1);
z1 = double(fixed.z) + shift(2);

m = 3; %Dimension of ODE
dt = 0.005;
tspan = 0:dt:100;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,m));

% Threshold value
eta = 0.05;

% Initial condition close to unstable orbit
x0 = [];
x0(1,:) = [x1+0.01; 0; z1];

% Controlled parameter
if abs(x0(1,1) - x1) + abs(x0(1,3) - z1) <= eta 
    mu(1) = mustar + C1(1)*(x0(1,1) - x1) + C1(2)*(x0(1,3) - z1);
else
    mu(1) = mustar;
end

% Initialize trajectory
[~,sol] = ode45(@(t,x) Sprott(x,mu(1)),tspan,x0(1,:));

% Initialize Controlled Solution
xc = [];
yc = [];
zc = [];

% Controlled orbit
kfinal = 50; % number of times the trajectory pierces Poincare section
for k = 2:kfinal
    
    for j = 1:length(sol(:,1))-1
       if  (sol(j,2) >= 0 && sol(j+1,2) < 0)  
            ind = j+1;
            
            % Controlled solution
            xc = [xc; sol(1:ind,1)];
            yc = [yc; sol(1:ind,2)];
            zc = [zc; sol(1:ind,3)];
            
            break
        end 
    end
   
    
    x0(k,:) = [sol(ind,1); sol(ind,2); sol(ind,3)];
    if abs(x0(k,1) - x1) + abs(x0(k,3) - z1) <= eta 
        mu(k) = mustar + C1(1)*(x0(k,1) - x1) + C1(2)*(x0(k,3) - z1);
    else
        mu(k) = mustar;
    end
    
    [~,sol] = ode45(@(t,x) Sprott(x,mu(k)),tspan,x0(k,:));
end

% Last Iteration of Controlled solution
xc = [xc; sol(1:ind,1)];
yc = [yc; sol(1:ind,2)];
zc = [zc; sol(1:ind,3)];

%% Identify period 2 orbit 

% Approximate (shifted) period 2 point
x21 = -0.6301;
z21 = 0.2206;
x22 = double(fs(x21,z21)); 
z22 = double(gs(x21,z21)); 

% Jacobian matrices
A21 = double([dfx(x21,z21,mu_s), dfz(x21,z21,mu_s); dgx(x21,z21,mu_s), dgz(x21,z21,mu_s)]);
B21 = double([dfp(x21,z21,mu_s); dgp(x21,z21,mu_s)]);
A22 = double([dfx(x22,z22,mu_s), dfz(x22,z22,mu_s); dgx(x22,z22,mu_s), dgz(x22,z22,mu_s)]);
B22 = double([dfp(x22,z22,mu_s); dgp(x22,z22,mu_s)]);

% Semidefinite constraints and optimization
M = [Q Q*A22' + Y'*B22'; A22*Q+B22*Y Q];
lmi = [M - I >= 0];
optimize(lmi)

% Control matrix for first period 2 point
C22 = value(Y)/(value(Q)); 

% First iterate matrix
A21 = (A22 + B22*C22)*A21;
B21 = (A22 + B22*C22)*B21;

% Semidefinite constraints and optimization
M = [Q Q*A21' + Y'*B21'; A21*Q+B21*Y Q];
lmi = [M - I >= 0];
optimize(lmi)

% Control matrix for second period 2 point
C21 = value(Y)/(value(Q));

%% Controlled trajectory (period 2 point)

% Shift variables back
x21 = -0.6301 + shift(1);
z21 = 0.2206 + shift(2);
x22 = double(fs(-0.6301,0.2206)) + shift(1);
z22 = double(gs(-0.6301,0.2206)) + shift(2);

% Threshold value
eta = 0.05;

% Initial condition close to unstable period 2 orbit
x0 = [];
x0(1,:) = [x21+0.01; 0; z21];

% Controlled parameter
if (x0(1,1) - x21)^2 + (x0(1,3) - z21)^2 <= eta 
    mu(1) = mustar + C21(1)*(x0(1,1) - x21) + C21(2)*(x0(1,3) - z21);
elseif (x0(1,1) - x22)^2 + (x0(1,3) - z22)^2 <= eta 
    mu(1) = mustar + C22(1)*(x0(1,1) - x22) + C22(2)*(x0(1,3) - z22);
else
    mu(1) = mustar;
end

% Initialize trajectory
[~,sol] = ode45(@(t,x) Sprott(x,mu(1)),tspan,x0(1,:));

% Initialize Controlled Solution
xc = [];
yc = [];
zc = [];

% Controlled orbit
kfinal = 50; % number of times the trajectory pierces Poincare section
for k = 2:kfinal
    
    for j = 1:length(sol(:,1))-1
       if  (sol(j,2) >= 0 && sol(j+1,2) < 0)  
            ind = j+1;
            
            % Controlled solution
            xc = [xc; sol(1:ind,1)];
            yc = [yc; sol(1:ind,2)];
            zc = [zc; sol(1:ind,3)];
            
            break
        end 
    end
   
    
    x0(k,:) = [sol(ind,1); sol(ind,2); sol(ind,3)];
    if (x0(k,1) - x21)^2 + (x0(k,3) - z21)^2 <= eta 
        mu(k) = mustar + C21(1)*(x0(k,1) - x21) + C21(2)*(x0(k,3) - z21);
    elseif (x0(k,1) - x22)^2 + (x0(k,3) - z22)^2 <= eta 
        mu(k) = mustar + C22(1)*(x0(k,1) - x22) + C22(2)*(x0(k,3) - z22);
    else
        mu(k) = mustar;
    end
    
    [~,sol] = ode45(@(t,x) Sprott(x,mu(k)),tspan,x0(k,:));
end

% Last Iteration of Controlled solution
xc = [xc; sol(1:ind,1)];
yc = [yc; sol(1:ind,2)];
zc = [zc; sol(1:ind,3)];

%% fsolve command to identify periodic orbits

per = 2;
x0 = rand(2,1)-0.5;

% option to display output and use Jacobian
options=optimset('Display','iter','Jacobian','on','MaxIter',10000,'Algorithm','levenberg-marquardt','TolFun',1e-6,'TolX',1e-6);

% call fsolve to find a periodic point of period = per above
xper = fsolve(@(x) Pmap(x,per,Xi_new,pow,mu_s),x0,options);

%% Sprott right-hand-side

function dx = Sprott(x,mu)
    
    dx = [x(2); x(3); -x(1) + x(2)^2 - mu*x(3)];

end

%% Discovered SINDy de Jong map 

function [F,J] = Pmap(x,per,Xi_new,pow,mu_s)

    % Initialize
    yold = x;

    % We are solving F(x) = x, or F(x) - x = 0, so the -x is the first part
    % of the map to put in
    F = -x;
    J = -eye(2);
    Anew = eye(2);

    % xnp1 = f(xn,zn) and znp1 = g(xn,zn)
    syms x z p
    f(x,z) = dot(prod([x z mu_s].^pow,2),Xi_new(:,1));
    g(x,z) = dot(prod([x z mu_s].^pow,2),Xi_new(:,2));

    % Symbolic derivatives
    dfx = diff(f,x);
    dfz = diff(f,z);
    dgx = diff(g,x);
    dgz = diff(g,z);
    
    for n = 1:per
       
        ynew(1) = double(f(yold(1),yold(2)));
        ynew(2) = double(g(yold(1),yold(2)));
        
        % Jacobian at each step
        Aold(1,1) = double(dfx(yold(1),yold(2))); 
        Aold(1,2) = double(dfz(yold(1),yold(2))); 
        Aold(2,1) = double(dgx(yold(1),yold(2))); 
        Aold(2,2) = double(dgz(yold(1),yold(2)));
        Anew = Aold*Anew;
        
        yold = ynew;
        
    end
    
    % Add in updates along the orbit
    F = F + ynew';
    J = J + Anew;
    
end





