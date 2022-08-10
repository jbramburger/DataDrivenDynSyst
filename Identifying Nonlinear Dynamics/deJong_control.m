% -------------------------------------------------------------------------
% Controlling Chaos in the de Jong Map
%
% We consider the two-dimensional de Jong map:
%
%     F1(x,y) = sin(a*y) - cos(b*x) 
%     F2(x,y) = sin(c*x) - cos(d*y) 
%
% with parameter values...
%
% TBD...
%
% This script accompanies Section 3.3 of Data-Driven Methods for
% Dynamic Systems. 
%
% This script uses YALMIP to translate the linear matrix inequality into a 
% semidefinite program to be solved by MATLAB. It can be downloaded at 
% https://yalmip.github.io/download/ and below it is assumed to be stored 
% in the folder 'YALMIP-master'. 
%
% Author: Jason J. Bramburger
% -------------------------------------------------------------------------

% Clean workspace
clear all; close all; clc

%% Add path to YALMIP (assumed to be in the folder 'YALMIP-master')

addpath(genpath('YALMIP-master'))

%% Initializations

% Map parameters
a = 0.970;
b = -1.899;
c = 1.381;
d = -1.506;

% Epsilon value
eps = 1e-3;
I = eps*eye(4);

% Number of iterations
N = 1e4;

% Threshold value (eta)
eta = 0.05;

% SDP variables
Q = sdpvar(2,2);
Y = sdpvar(1,2,'full');

%% Fixed point stabilization

% Use symbolic variables to identify the unique fixed point
syms x y
f(x,y) = sin(a*y) - cos(b*x);
g(x,y) = sin(c*x) - cos(d*y);
sol1 = vpasolve([f(x,y) == x, g(x,y) == y], [x, y]);
x1 = double(sol1.x);
y1 = double(sol1.y);

% Linearization matrices
A = [b*sin(b*x1), a*cos(a*y1); c*cos(c*x1), d*sin(d*y1)];
B = [y1*cos(a*y1); 0];

% Semidefinite constraints and optimization
M = [Q Q*A' + Y'*B'; A*Q+B*Y Q];
lmi = [M - I >= 0];
optimize(lmi)

% Control matrix for fixed points
K1 = value(Y)/(value(Q)); 

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
        a_control = a + K1*[xc(n) - x1; yc(n) - y1]; % parameter change
        
        xc(n+1) = sin(a_control*yc(n)) - cos(b*xc(n));
        yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
    % Otherwise continue with de Jong iterations
    else
        xc(n+1) = sin(a*yc(n)) - cos(b*xc(n));
        yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
    end
end

%Plot controlled x trajectory
figure(1)
plot(xc(1:1000),'.-','Color',[1 69/255 79/255],'MarkerSize',10)
set(gca,'FontSize',16)

%Plot controlled y trajectory
figure(2)
plot(yc(1:1000),'.-','Color',[36/255 122/255 254/255],'MarkerSize',10)
set(gca,'FontSize',16)

%% Period 2 stabilization

% Use symbolic variables to identify the a point on the period 2 orbit
syms x y
f(x,y) = sin(a*y) - cos(b*x);
g(x,y) = sin(c*x) - cos(d*y);
sol2 = vpasolve([f(f(x,y),g(x,y)) == x, g(f(x,y),g(x,y)) == y], [x, y]);
x21 = double(sol2.x);
y21 = double(sol2.y);
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
K21 = value(Y)/(value(Q)); 

% Second iterate matrix
A22 = (A21 + B21*K21)*A22;
B22 = (A21 + B21*K21)*B22;

% Semidefinite constraints and optimization
M = [Q Q*A22' + Y'*B22'; A22*Q+B22*Y Q];
lmi = [M - I >= 0];
optimize(lmi)

% Control matrix for second period 2 point
K22 = value(Y)/(value(Q));

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
        a_control = a + K21*[xc(n) - x21; yc(n) - y21]; % parameter change
        
        xc(n+1) = sin(a_control*yc(n)) - cos(b*xc(n));
        yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
        
    elseif (xc(n) - x22)^2 + (yc(n) - y22)^2 <= eta %Kick system when close to other period 2 point
        a_control = a + K22*[xc(n) - x22; yc(n) - y22]; % parameter change
        
        xc(n+1) = sin(a_control*yc(n)) - cos(b*xc(n));
        yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
    % Otherwise continue with de Jong iterations
    else
        xc(n+1) = sin(a*yc(n)) - cos(b*xc(n));
        yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
    end
end

%Plot controlled x trajectory
figure(3)
plot(xc(1:1000),'.-','Color',[1 69/255 79/255],'MarkerSize',10)
set(gca,'FontSize',16)

%Plot controlled y trajectory
figure(4)
plot(yc(1:1000),'.-','Color',[36/255 122/255 254/255],'MarkerSize',10)
set(gca,'FontSize',16)

%% Period 6 stabilization

% Period 6 point is identified with code at the end of the script
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
K71 = value(Y)/(value(Q)); 

% Second iterate matrix
A72 = (A71 + B71*K71)*A72;
B72 = (A71 + B71*K71)*B72;

% Semidefinite constraints and optimization
M = [Q Q*A72' + Y'*B72'; A72*Q+B72*Y Q];
lmi = [M - I >= 0];
optimize(lmi)

% Control matrix for second period 7 point
K72 = value(Y)/(value(Q));

% Third iterate matrix
A73 = (A72 + B72*K72)*A73;
B73 = (A72 + B72*K72)*B73;

% Semidefinite constraints and optimization
M = [Q Q*A73' + Y'*B73'; A73*Q+B73*Y Q];
lmi = [M - I >= 0];
optimize(lmi)

% Control matrix for third period 7 point
K73 = value(Y)/(value(Q));

% Fourth iterate matrix
A74 = (A73 + B73*K73)*A74;
B74 = (A73 + B73*K73)*B74;

% Semidefinite constraints and optimization
M = [Q Q*A74' + Y'*B74'; A74*Q+B74*Y Q];
lmi = [M - I >= 0];
optimize(lmi)

% Control matrix for fourth period 7 point
K74 = value(Y)/(value(Q));

% Fifth iterate matrix
A75 = (A74 + B74*K74)*A75;
B75 = (A74 + B74*K74)*B75;

% Semidefinite constraints and optimization
M = [Q Q*A75' + Y'*B75'; A75*Q+B75*Y Q];
lmi = [M - I >= 0];
optimize(lmi)

% Control matrix for fifth period 7 point
K75 = value(Y)/(value(Q));

% Sixth iterate matrix
A76 = (A75 + B75*K75)*A76;
B76 = (A75 + B75*K75)*B76;

% Semidefinite constraints and optimization
M = [Q Q*A76' + Y'*B76'; A76*Q+B76*Y Q];
lmi = [M - I >= 0];
optimize(lmi)

% Control matrix for sixth period 7 point
K76 = value(Y)/(value(Q));

% Seventh iterate matrix
A77 = (A76 + B76*K76)*A77;
B77 = (A76 + B76*K76)*B77;

% Semidefinite constraints and optimization
M = [Q Q*A77' + Y'*B77'; A77*Q+B77*Y Q];
lmi = [M - I >= 0];
optimize(lmi)

% Control matrix for seventh period 7 point
K77 = value(Y)/(value(Q));

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
        a_control = a + K71*[xc(n) - x71; yc(n) - y71]; % parameter change
        
        xc(n+1) = sin(a_control*yc(n)) - cos(b*xc(n));
        yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
        
    elseif (xc(n) - x72)^2 + (yc(n) - y72)^2 <= eta %Kick system when close to second period 3 point
        a_control = a + K72*[xc(n) - x72; yc(n) - y72]; % parameter change
        
        xc(n+1) = sin(a_control*yc(n)) - cos(b*xc(n));
        yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
    
    elseif (xc(n) - x73)^2 + (yc(n) - y73)^2 <= eta %Kick system when close to second period 3 point
        a_control = a + K73*[xc(n) - x73; yc(n) - y73]; % parameter change
        
        xc(n+1) = sin(a_control*yc(n)) - cos(b*xc(n));
        yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
        
    elseif (xc(n) - x74)^2 + (yc(n) - y74)^2 <= eta %Kick system when close to second period 3 point
        a_control = a + K74*[xc(n) - x74; yc(n) - y74]; % parameter change
        
        xc(n+1) = sin(a_control*yc(n)) - cos(b*xc(n));
        yc(n+1) = sin(c*xc(n)) - cos(d*yc(n)); 
    
    elseif (xc(n) - x75)^2 + (yc(n) - y75)^2 <= eta %Kick system when close to second period 3 point
        a_control = a + K75*[xc(n) - x75; yc(n) - y75]; % parameter change
        
        xc(n+1) = sin(a_control*yc(n)) - cos(b*xc(n));
        yc(n+1) = sin(c*xc(n)) - cos(d*yc(n)); 
        
    elseif (xc(n) - x76)^2 + (yc(n) - y76)^2 <= eta %Kick system when close to second period 3 point
        a_control = a + K76*[xc(n) - x76; yc(n) - y76]; % parameter change
        
        xc(n+1) = sin(a_control*yc(n)) - cos(b*xc(n));
        yc(n+1) = sin(c*xc(n)) - cos(d*yc(n)); 
        
    elseif (xc(n) - x77)^2 + (yc(n) - y77)^2 <= eta %Kick system when close to second period 3 point
        a_control = a + K77*[xc(n) - x77; yc(n) - y77]; % parameter change
        
        xc(n+1) = sin(a_control*yc(n)) - cos(b*xc(n));
        yc(n+1) = sin(c*xc(n)) - cos(d*yc(n)); 
    % Otherwise continue with de Jong iterations
    else
        xc(n+1) = sin(a*yc(n)) - cos(b*xc(n));
        yc(n+1) = sin(c*xc(n)) - cos(d*yc(n));
    end
end

%Plot controlled x trajectory
figure(5)
plot(xc(1:1000),'.-','Color',[1 69/255 79/255],'MarkerSize',10)
set(gca,'FontSize',16)

%Plot controlled y trajectory
figure(6)
plot(yc(1:1000),'.-','Color',[36/255 122/255 254/255],'MarkerSize',10)
set(gca,'FontSize',16)

%% The following can be used to identify periodic orbits of the map

per = 7;
x0 = 2*rand(2,1)-1;

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


