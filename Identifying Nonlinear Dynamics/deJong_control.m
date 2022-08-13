% -------------------------------------------------------------------------
% Stabilizing Periodic Orbits in the de Jong Map
%
% In this script we implement the method of linear matrix inequalities
% (LMIs) to stabilize periodic orbits of the two-dimensional de Jong map:
%
%     F1(x,y) = sin(a*y) - cos(b*x) 
%     F2(x,y) = sin(c*x) - cos(d*y) 
%
% with parameter values (b,c,d) = (-1.899,1.381,-1.506). The control 
% parameter is a, centred around the value a = 0.970 at which the de Jong 
% map is chaotic. 
%
% Below we stabilize a fixed point, period 2 orbit, and a period 7 orbit.
% We also include a switching controller that stabilizes the period 7
% orbit, then the fixed point, and then the period 2 orbit. At the end of
% this script is a root-finding procedure to identify periodic points of
% the de Jong map. This procedure uses the built-in MATLAB function fsolve,
% which is included with the optimization toolbox. This root-finding
% procedure is not crucial to running the control processes. However, we do
% use symbolic variables to identify the fixed and period 2 points of the 
% map, requiring the symbolic math toolbox for MATLAB. The resulting values 
% from the symbolic computation is provided in the comments below in case 
% one does not have the symbolic math toolbox. 
%
% This script accompanies Section 3.3 of Data-Driven Methods for
% Dynamic Systems. 
%
% This script uses YALMIP to translate the linear matrix inequality into a 
% semidefinite program to be solved by a variety of different solvers. 
% This work was undertaken originally with the solver MOSEK, although most
% others (such as SEDUMI or SDPT3) will work as well. For a complete 
% discussion of solvers to use with YALMIP please visit: 
%      https://yalmip.github.io/allsolvers/
% 
% Both YALMIP and MOSEK are freely available and can be downloaded at:
%     YALMIP: https://yalmip.github.io/download/ 
%     MOSEK: https://www.mosek.com/
%
% Author: Jason J. Bramburger
% -------------------------------------------------------------------------

% Clean workspace
clear all; close all; clc

%% Add path to YALMIP (assumed to be in the folder 'YALMIP-master')

addpath(genpath('YALMIP-master'))

%% Initializations

% Figure counter
figCount = 1;

% Map parameters
a = 0.970;
b = -1.899;
c = 1.381;
d = -1.506;

% Epsilon value for numerical strict inequalities
eps = 1e-3;
I = eps*eye(4);

% Number of iterations
N = 1e4;

% Threshold value (eta)
eta = 0.05; % <---- this is eta^2 in the text

% SDP variables
Q = sdpvar(2,2);
Y = sdpvar(1,2,'full');

%% Fixed point stabilization

% Use symbolic variables to identify the unique fixed point
%   (values of x1 and y1 are provided below if symbolic math toolbox is
%   not installed)
syms x y
f(x,y) = sin(a*y) - cos(b*x);
g(x,y) = sin(c*x) - cos(d*y);
sol1 = vpasolve([f(x,y) == x, g(x,y) == y], [x, y]);
x1 = double(sol1.x); % value is -0.8243
y1 = double(sol1.y); % value is -0.9891

% Linearization matrices
A = [b*sin(b*x1), a*cos(a*y1); c*cos(c*x1), d*sin(d*y1)];
B = [y1*cos(a*y1); 0];

% Semidefinite constraints and optimization
M = [Q Q*A' + Y'*B'; A*Q+B*Y Q];
lmi = [M - I >= 0];
optimize(lmi)

% Control matrix for fixed points
C1 = value(Y)/(value(Q)); 

%% Control output to a fixed point

% Controlled variables
xc = zeros(N+1,1);
yc = zeros(N+1,1);
xc(1) = 0.1;
yc(1) = 0.0;

% Controlled iterations
for n = 1:N
    
    % If close, then control
    if (xc(n) - x1)^2 + (yc(n) - y1)^2 <= eta %Kick system when close to fixed points
        a_control = a + C1*[xc(n) - x1; yc(n) - y1]; % parameter change
        
        xc(n+1) = sin(a_control*yc(n)) - cos(b*xc(n));
        yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
    % Otherwise continue with de Jong iterations
    else
        xc(n+1) = sin(a*yc(n)) - cos(b*xc(n));
        yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
    end
end

%Plot controlled x trajectory
figure(figCount)
plot(xc(1:1000),'.-','Color',[1 69/255 79/255],'MarkerSize',10)
set(gca,'FontSize',16)
figCount = figCount + 1;

%Plot controlled y trajectory
figure(figCount)
plot(yc(1:1000),'.-','Color',[36/255 122/255 254/255],'MarkerSize',10)
set(gca,'FontSize',16)
figCount = figCount + 1;

%% Period 2 stabilization

% Use symbolic variables to identify the a point on the period 2 orbit
%   (values of x21 and y21 are provided below if symbolic math toolbox is
%   not installed)
syms x y
f(x,y) = sin(a*y) - cos(b*x);
g(x,y) = sin(c*x) - cos(d*y);
sol2 = vpasolve([f(f(x,y),g(x,y)) == x, g(f(x,y),g(x,y)) == y], [x, y]);
x21 = double(sol2.x); % value is -0.5515 
y21 = double(sol2.y); % value is -0.5963
x22 = double(f(x21,y21));
y22 = double(g(x21,y21));

% Linearization matrices
A21 = [b*sin(b*x21), a*cos(a*y21); c*cos(c*x21), d*sin(d*y21)];
B21 = [y21*cos(a*y21); 0];
A22 = [b*sin(b*x22), a*cos(a*y22); c*cos(c*x22), d*sin(d*y22)];
B22 = [y22*cos(a*y22); 0];

% Semidefinite constraints and optimization
M = [Q Q*A21' + Y'*B21'; A21*Q+B21*Y Q];
lmi = [M - I >= 0];
optimize(lmi)

% Control matrix for first period 2 point
C21 = value(Y)/(value(Q)); 

% Second iterate matrix
A22 = (A21 + B21*C21)*A22;
B22 = (A21 + B21*C21)*B22;

% Semidefinite constraints and optimization
M = [Q Q*A22' + Y'*B22'; A22*Q+B22*Y Q];
lmi = [M - I >= 0];
optimize(lmi)

% Control matrix for second period 2 point
C22 = value(Y)/(value(Q));

%% Control output to a period 2 orbit

% Controlled variables
xc = zeros(N+1,1);
yc = zeros(N+1,1);
xc(1) = 0.1;
yc(1) = 0.0;

% Controlled iterations
for n = 1:N
    
    % If close, then control
    if (xc(n) - x21)^2 + (yc(n) - y21)^2 <= eta %Kick system when close to first period 2 point
        a_control = a + C21*[xc(n) - x21; yc(n) - y21]; % parameter change
        
        xc(n+1) = sin(a_control*yc(n)) - cos(b*xc(n));
        yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
        
    elseif (xc(n) - x22)^2 + (yc(n) - y22)^2 <= eta %Kick system when close to other period 2 point
        a_control = a + C22*[xc(n) - x22; yc(n) - y22]; % parameter change
        
        xc(n+1) = sin(a_control*yc(n)) - cos(b*xc(n));
        yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
    % Otherwise continue with de Jong iterations
    else
        xc(n+1) = sin(a*yc(n)) - cos(b*xc(n));
        yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
    end
end

%Plot controlled x trajectory
figure(figCount)
plot(xc(1:1000),'.-','Color',[1 69/255 79/255],'MarkerSize',10)
set(gca,'FontSize',16)
figCount = figCount + 1;

%Plot controlled y trajectory
figure(figCount)
plot(yc(1:1000),'.-','Color',[36/255 122/255 254/255],'MarkerSize',10)
set(gca,'FontSize',16)
figCount = figCount + 1;

%% Period 7 stabilization

% Period 7 point is identified with code at the end of the script
x71 = 1.0219;
y71 = 0.0382;
x72 = double(f(x71,y71));
y72 = double(g(x71,y71));
x73 = double(f(x72,y72));
y73 = double(g(x72,y72));
x74 = double(f(x73,y73));
y74 = double(g(x73,y73));
x75 = double(f(x74,y74));
y75 = double(g(x74,y74));
x76 = double(f(x75,y75));
y76 = double(g(x75,y75));
x77 = double(f(x76,y76));
y77 = double(g(x76,y76));

% Linearization matrices
A71 = [b*sin(b*x71), a*cos(a*y71); c*cos(c*x71), d*sin(d*y71)];
B71 = [y71*cos(a*y71); 0];
A72 = [b*sin(b*x72), a*cos(a*y72); c*cos(c*x72), d*sin(d*y72)];
B72 = [y72*cos(a*y72); 0];
A73 = [b*sin(b*x73), a*cos(a*y73); c*cos(c*x73), d*sin(d*y73)];
B73 = [y73*cos(a*y73); 0];
A74 = [b*sin(b*x74), a*cos(a*y74); c*cos(c*x74), d*sin(d*y74)];
B74 = [y74*cos(a*y74); 0];
A75 = [b*sin(b*x75), a*cos(a*y75); c*cos(c*x75), d*sin(d*y75)];
B75 = [y75*cos(a*y75); 0];
A76 = [b*sin(b*x76), a*cos(a*y76); c*cos(c*x76), d*sin(d*y76)];
B76 = [y76*cos(a*y76); 0];
A77 = [b*sin(b*x77), a*cos(a*y77); c*cos(c*x77), d*sin(d*y77)];
B77 = [y77*cos(a*y77); 0];

% Semidefinite constraints and optimization
M = [Q Q*A71' + Y'*B71'; A71*Q+B71*Y Q];
lmi = [M - I >= 0];
optimize(lmi)

% Control matrix for first period 7 point
C71 = value(Y)/(value(Q)); 

% Second iterate matrix
A72 = (A71 + B71*C71)*A72;
B72 = (A71 + B71*C71)*B72;

% Semidefinite constraints and optimization
M = [Q Q*A72' + Y'*B72'; A72*Q+B72*Y Q];
lmi = [M - I >= 0];
optimize(lmi)

% Control matrix for second period 7 point
C72 = value(Y)/(value(Q));

% Third iterate matrix
A73 = (A72 + B72*C72)*A73;
B73 = (A72 + B72*C72)*B73;

% Semidefinite constraints and optimization
M = [Q Q*A73' + Y'*B73'; A73*Q+B73*Y Q];
lmi = [M - I >= 0];
optimize(lmi)

% Control matrix for third period 7 point
C73 = value(Y)/(value(Q));

% Fourth iterate matrix
A74 = (A73 + B73*C73)*A74;
B74 = (A73 + B73*C73)*B74;

% Semidefinite constraints and optimization
M = [Q Q*A74' + Y'*B74'; A74*Q+B74*Y Q];
lmi = [M - I >= 0];
optimize(lmi)

% Control matrix for fourth period 7 point
C74 = value(Y)/(value(Q));

% Fifth iterate matrix
A75 = (A74 + B74*C74)*A75;
B75 = (A74 + B74*C74)*B75;

% Semidefinite constraints and optimization
M = [Q Q*A75' + Y'*B75'; A75*Q+B75*Y Q];
lmi = [M - I >= 0];
optimize(lmi)

% Control matrix for fifth period 7 point
C75 = value(Y)/(value(Q));

% Sixth iterate matrix
A76 = (A75 + B75*C75)*A76;
B76 = (A75 + B75*C75)*B76;

% Semidefinite constraints and optimization
M = [Q Q*A76' + Y'*B76'; A76*Q+B76*Y Q];
lmi = [M - I >= 0];
optimize(lmi)

% Control matrix for sixth period 7 point
C76 = value(Y)/(value(Q));

% Seventh iterate matrix
A77 = (A76 + B76*C76)*A77;
B77 = (A76 + B76*C76)*B77;

% Semidefinite constraints and optimization
M = [Q Q*A77' + Y'*B77'; A77*Q+B77*Y Q];
lmi = [M - I >= 0];
optimize(lmi)

% Control matrix for seventh period 7 point
C77 = value(Y)/(value(Q));

%% Control output to a period 7 orbit

% Controlled variables
xc = zeros(N+1,1);
yc = zeros(N+1,1);
xc(1) = 0.1;
yc(1) = 0.0;

% Controlled iterations
for n = 1:N
    
    % If close, then control
    if (xc(n) - x71)^2 + (yc(n) - y71)^2 <= eta %Kick system when close to first period 3 point
        a_control = a + C71*[xc(n) - x71; yc(n) - y71]; % parameter change
        
        xc(n+1) = sin(a_control*yc(n)) - cos(b*xc(n));
        yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
        
    elseif (xc(n) - x72)^2 + (yc(n) - y72)^2 <= eta %Kick system when close to second period 3 point
        a_control = a + C72*[xc(n) - x72; yc(n) - y72]; % parameter change
        
        xc(n+1) = sin(a_control*yc(n)) - cos(b*xc(n));
        yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
    
    elseif (xc(n) - x73)^2 + (yc(n) - y73)^2 <= eta %Kick system when close to second period 3 point
        a_control = a + C73*[xc(n) - x73; yc(n) - y73]; % parameter change
        
        xc(n+1) = sin(a_control*yc(n)) - cos(b*xc(n));
        yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
        
    elseif (xc(n) - x74)^2 + (yc(n) - y74)^2 <= eta %Kick system when close to second period 3 point
        a_control = a + C74*[xc(n) - x74; yc(n) - y74]; % parameter change
        
        xc(n+1) = sin(a_control*yc(n)) - cos(b*xc(n));
        yc(n+1) = sin(c*xc(n)) - cos(d*yc(n)); 
    
    elseif (xc(n) - x75)^2 + (yc(n) - y75)^2 <= eta %Kick system when close to second period 3 point
        a_control = a + C75*[xc(n) - x75; yc(n) - y75]; % parameter change
        
        xc(n+1) = sin(a_control*yc(n)) - cos(b*xc(n));
        yc(n+1) = sin(c*xc(n)) - cos(d*yc(n)); 
        
    elseif (xc(n) - x76)^2 + (yc(n) - y76)^2 <= eta %Kick system when close to second period 3 point
        a_control = a + C76*[xc(n) - x76; yc(n) - y76]; % parameter change
        
        xc(n+1) = sin(a_control*yc(n)) - cos(b*xc(n));
        yc(n+1) = sin(c*xc(n)) - cos(d*yc(n)); 
        
    elseif (xc(n) - x77)^2 + (yc(n) - y77)^2 <= eta %Kick system when close to second period 3 point
        a_control = a + C77*[xc(n) - x77; yc(n) - y77]; % parameter change
        
        xc(n+1) = sin(a_control*yc(n)) - cos(b*xc(n));
        yc(n+1) = sin(c*xc(n)) - cos(d*yc(n)); 
    % Otherwise continue with de Jong iterations
    else
        xc(n+1) = sin(a*yc(n)) - cos(b*xc(n));
        yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
    end
end

%Plot controlled x trajectory
figure(figCount)
plot(xc(1:1000),'.-','Color',[1 69/255 79/255],'MarkerSize',10)
set(gca,'FontSize',16)
figCount = figCount + 1;

%Plot controlled y trajectory
figure(figCount)
plot(yc(1:1000),'.-','Color',[36/255 122/255 254/255],'MarkerSize',10)
set(gca,'FontSize',16)
figCount = figCount + 1;

%% Plotting the attractor and embedded periodic points

% Period 8 points
x81 = -0.9756;
y81 = -1.0447;
x82 = double(f(x81,y81));
y82 = double(g(x81,y81));
x83 = double(f(x82,y82));
y83 = double(g(x82,y82));
x84 = double(f(x83,y83));
y84 = double(g(x83,y83));
x85 = double(f(x84,y84));
y85 = double(g(x84,y84));
x86 = double(f(x85,y85));
y86 = double(g(x85,y85));
x87 = double(f(x86,y86));
y87 = double(g(x86,y86));
x88 = double(f(x87,y87));
y88 = double(g(x87,y87));

% Period 9 points
x91 = 0.7520;
y91 = -0.0174;
x92 = double(f(x91,y91));
y92 = double(g(x91,y91));
x93 = double(f(x92,y92));
y93 = double(g(x92,y92));
x94 = double(f(x93,y93));
y94 = double(g(x93,y93));
x95 = double(f(x94,y94));
y95 = double(g(x94,y94));
x96 = double(f(x95,y95));
y96 = double(g(x95,y95));
x97 = double(f(x96,y96));
y97 = double(g(x96,y96));
x98 = double(f(x97,y97));
y98 = double(g(x97,y97));
x99 = double(f(x98,y98));
y99 = double(g(x98,y98));

% Period 10 points
x101 = 0.0562;
y101 = 0.2661;
x102 = double(f(x101,y101));
y102 = double(g(x101,y101));
x103 = double(f(x102,y102));
y103 = double(g(x102,y102));
x104 = double(f(x103,y103));
y104 = double(g(x103,y103));
x105 = double(f(x104,y104));
y105 = double(g(x104,y104));
x106 = double(f(x105,y105));
y106 = double(g(x105,y105));
x107 = double(f(x106,y106));
y107 = double(g(x106,y106));
x108 = double(f(x107,y107));
y108 = double(g(x107,y107));
x109 = double(f(x108,y108));
y109 = double(g(x108,y108));
x1010 = double(f(x109,y109));
y1010 = double(g(x109,y109));

% Initial attractor points (IC = (0,0))
xplot = zeros(1e5,1);
yplot = zeros(1e5,1);

% Generate trajectory data
for n = 2:1e5
    xplot(n) = sin(a*yplot(n-1)) - cos(b*xplot(n-1));
    yplot(n) = sin(c*xplot(n-1)) - cos(d*yplot(n-1));
end

% Plot attractor
figure(figCount)
plot(xplot(100:end),yplot(100:end),'k.','MarkerSize',2)
hold on
plot(x1,y1,'.','Color',[1 69/255 79/255],'MarkerSize',30) % Fixed point

plot(x21,y21,'.','Color',[36/255 122/255 254/255],'MarkerSize',30) % Period 2 points
plot(x22,y22,'.','Color',[36/255 122/255 254/255],'MarkerSize',30)

plot(x71,y71,'.','Color',[146/255 208/255 80/255],'MarkerSize',30) % Period 7 points
plot(x72,y72,'.','Color',[146/255 208/255 80/255],'MarkerSize',30)
plot(x73,y73,'.','Color',[146/255 208/255 80/255],'MarkerSize',30)
plot(x74,y74,'.','Color',[146/255 208/255 80/255],'MarkerSize',30)
plot(x75,y75,'.','Color',[146/255 208/255 80/255],'MarkerSize',30)
plot(x76,y76,'.','Color',[146/255 208/255 80/255],'MarkerSize',30)
plot(x77,y77,'.','Color',[146/255 208/255 80/255],'MarkerSize',30)

plot(x81,y81,'.','Color',[254/255 174/255 0/255],'MarkerSize',30) % Period 8 poitns
plot(x82,y82,'.','Color',[254/255 174/255 0/255],'MarkerSize',30)
plot(x83,y83,'.','Color',[254/255 174/255 0/255],'MarkerSize',30)
plot(x84,y84,'.','Color',[254/255 174/255 0/255],'MarkerSize',30)
plot(x85,y85,'.','Color',[254/255 174/255 0/255],'MarkerSize',30)
plot(x86,y86,'.','Color',[254/255 174/255 0/255],'MarkerSize',30)
plot(x87,y87,'.','Color',[254/255 174/255 0/255],'MarkerSize',30)
plot(x88,y88,'.','Color',[254/255 174/255 0/255],'MarkerSize',30)

plot(x91,y91,'.','Color',[255/255 66/255 161/255],'MarkerSize',30) % Period 9 points
plot(x92,y92,'.','Color',[255/255 66/255 161/255],'MarkerSize',30)
plot(x93,y93,'.','Color',[255/255 66/255 161/255],'MarkerSize',30)
plot(x94,y94,'.','Color',[255/255 66/255 161/255],'MarkerSize',30)
plot(x95,y95,'.','Color',[255/255 66/255 161/255],'MarkerSize',30)
plot(x96,y96,'.','Color',[255/255 66/255 161/255],'MarkerSize',30)
plot(x97,y97,'.','Color',[255/255 66/255 161/255],'MarkerSize',30)
plot(x98,y98,'.','Color',[255/255 66/255 161/255],'MarkerSize',30)
plot(x99,y99,'.','Color',[255/255 66/255 161/255],'MarkerSize',30)

plot(x101,y101,'.','Color',[22/255 231/255 207/255],'MarkerSize',30) % Period 10 points
plot(x102,y102,'.','Color',[22/255 231/255 207/255],'MarkerSize',30)
plot(x103,y103,'.','Color',[22/255 231/255 207/255],'MarkerSize',30)
plot(x104,y104,'.','Color',[22/255 231/255 207/255],'MarkerSize',30)
plot(x105,y105,'.','Color',[22/255 231/255 207/255],'MarkerSize',30)
plot(x106,y106,'.','Color',[22/255 231/255 207/255],'MarkerSize',30)
plot(x107,y107,'.','Color',[22/255 231/255 207/255],'MarkerSize',30)
plot(x108,y108,'.','Color',[22/255 231/255 207/255],'MarkerSize',30)
plot(x109,y109,'.','Color',[22/255 231/255 207/255],'MarkerSize',30)
plot(x1010,y1010,'.','Color',[22/255 231/255 207/255],'MarkerSize',30)

xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
set(gca,'FontSize',16,'Xlim',[-2.1,2.1],'Ylim',[-2.1,2.1])
figCount = figCount + 1;

%% Switching between controlled orbits (7-->1-->2)

% Controlled variables
xc = zeros(N+1,1);
yc = zeros(N+1,1);
xc(1) = 0.5;
yc(1) = 0.0;

max1 = 250; % maximal number of iterates to stay near period 1 orbit
max7 = 250; % maximal number of iterates to stay near period 7 orbit
perOff = 100; % provides separation between 1 --> 2 transition for visualization
per1count = 1;
per7count = 1;
perOffcount = 1; 

% Controlled iterations
for n = 1:N
    
    if per7count <= max7
        % If close, then control
        if (xc(n) - x71)^2 + (yc(n) - y71)^2 <= eta %Kick system when close to first period 3 point
            a_control = a + C71*[xc(n) - x71; yc(n) - y71]; % parameter change

            xc(n+1) = sin(a_control*yc(n)) - cos(b*xc(n));
            yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
            per7count = per7count + 1;
            
        elseif (xc(n) - x72)^2 + (yc(n) - y72)^2 <= eta %Kick system when close to second period 3 point
            a_control = a + C72*[xc(n) - x72; yc(n) - y72]; % parameter change

            xc(n+1) = sin(a_control*yc(n)) - cos(b*xc(n));
            yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
            per7count = per7count + 1;

        elseif (xc(n) - x73)^2 + (yc(n) - y73)^2 <= eta %Kick system when close to second period 3 point
            a_control = a + C73*[xc(n) - x73; yc(n) - y73]; % parameter change

            xc(n+1) = sin(a_control*yc(n)) - cos(b*xc(n));
            yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
            per7count = per7count + 1;

        elseif (xc(n) - x74)^2 + (yc(n) - y74)^2 <= eta %Kick system when close to second period 3 point
            a_control = a + C74*[xc(n) - x74; yc(n) - y74]; % parameter change

            xc(n+1) = sin(a_control*yc(n)) - cos(b*xc(n));
            yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
            per7count = per7count + 1;

        elseif (xc(n) - x75)^2 + (yc(n) - y75)^2 <= eta %Kick system when close to second period 3 point
            a_control = a + C75*[xc(n) - x75; yc(n) - y75]; % parameter change

            xc(n+1) = sin(a_control*yc(n)) - cos(b*xc(n));
            yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
            per7count = per7count + 1;

        elseif (xc(n) - x76)^2 + (yc(n) - y76)^2 <= eta %Kick system when close to second period 3 point
            a_control = a + C76*[xc(n) - x76; yc(n) - y76]; % parameter change

            xc(n+1) = sin(a_control*yc(n)) - cos(b*xc(n));
            yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
            per7count = per7count + 1;

        elseif (xc(n) - x77)^2 + (yc(n) - y77)^2 <= eta %Kick system when close to second period 3 point
            a_control = a + C77*[xc(n) - x77; yc(n) - y77]; % parameter change

            xc(n+1) = sin(a_control*yc(n)) - cos(b*xc(n));
            yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
            per7count = per7count + 1;
        % Otherwise continue with de Jong iterations
        else
            xc(n+1) = sin(a*yc(n)) - cos(b*xc(n));
            yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
        end
    elseif per7count > max7 && per1count <= max1    
        % If close, then control
        if (xc(n) - x1)^2 + (yc(n) - y1)^2 <= eta %Kick system when close to fixed points
            a_control = a + C1*[xc(n) - x1; yc(n) - y1]; % parameter change

            xc(n+1) = sin(a_control*yc(n)) - cos(b*xc(n));
            yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
            per1count = per1count + 1;
        % Otherwise continue with de Jong iterations
        else
            xc(n+1) = sin(a*yc(n)) - cos(b*xc(n));
            yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
        end
    elseif per7count > max7 && per1count > max1 && perOffcount <= perOff
        xc(n+1) = sin(a*yc(n)) - cos(b*xc(n));
        yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
        perOffcount = perOffcount + 1;
        
    elseif per7count > max7 && per1count > max1 && perOffcount > perOff
        % If close, then control
        if (xc(n) - x21)^2 + (yc(n) - y21)^2 <= eta %Kick system when close to first period 2 point
            a_control = a + C21*[xc(n) - x21; yc(n) - y21]; % parameter change

            xc(n+1) = sin(a_control*yc(n)) - cos(b*xc(n));
            yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));

        elseif (xc(n) - x22)^2 + (yc(n) - y22)^2 <= eta %Kick system when close to other period 2 point
            a_control = a + C22*[xc(n) - x22; yc(n) - y22]; % parameter change

            xc(n+1) = sin(a_control*yc(n)) - cos(b*xc(n));
            yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
        % Otherwise continue with de Jong iterations
        else
            xc(n+1) = sin(a*yc(n)) - cos(b*xc(n));
            yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
        end
    end
end

%Plot controlled x trajectory
figure(figCount)
plot(xc(1:1200),'.-','Color',[1 69/255 79/255],'MarkerSize',10)
set(gca,'FontSize',16)
figCount = figCount + 1;

%Plot controlled y trajectory
figure(figCount)
plot(yc(1:1200),'.-','Color',[36/255 122/255 254/255],'MarkerSize',10)
set(gca,'FontSize',16)
figCount = figCount + 1;


%% The following can be used to identify periodic orbits of the map

per = 9;
x0 = 4*rand(2,1)-2;

% option to display output and use Jacobian
options=optimset('Display','iter','Jacobian','on','MaxIter',10000,'Algorithm','levenberg-marquardt','TolFun',1e-15,'TolX',1e-15);

% call fsolve to find a periodic point of period = per above
xper = fsolve(@(x) map(x,per),x0,options);

%% Iterated de Jong map 

function [F,J] = map(x,per)

    % de Jong map parameters
    a = 0.970;
    b = -1.899;
    c = 1.381;
    d = -1.506;

    % Initialize
    yold = x;

    % We are solving F(x) = x, or F(x) - x = 0, so the -x is the first part
    % of the map to put in
    F = -x;
    J = -eye(2);
    Anew = eye(2);
    
    for n = 1:per
       
        % Map iterates
        ynew(1) = sin(a*yold(2)) - cos(b*yold(1));
        ynew(2) = sin(c*yold(1)) - cos(d*yold(2));
        
        % Jacobian at each step
        Aold(1,1) = b*sin(b*yold(1)); 
        Aold(1,2) = a*cos(a*yold(2)); 
        Aold(2,1) = c*cos(c*yold(1)); 
        Aold(2,2) = d*sin(d*yold(2));
        Anew = Aold*Anew;
        
        yold = ynew;
        
    end
    
    % Add in updates along the orbit
    F = F + ynew';
    J = J + Anew;
    
end


