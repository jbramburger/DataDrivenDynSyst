% -------------------------------------------------------------------------
% Delay Coordinate Dynamic Mode Decomposition 
%
% This code applies the delay-coordinate DMD approach to scalar time series.
% There are two examples included: one synthetic example generated from a
% simple linear span of sines and cosines with differing frequencies and 
% another from the measurements of the nonlinear oscillations of the Van 
% der Pol oscillator. We further take the SVD of the Hankel matrix 
% associated to the Van der Pol measurements to demonstrate the nearly 
% linear oscillation of the observables coming from the left singular
% vectors in U.
%
% This script accompanies Section 2.1.1 of Computational Methods for
% Dynamical Systems. 
%
% Author: Jason J. Bramburger
% -------------------------------------------------------------------------

% Clean workspace
clear all; close all; clc

%% Simple synthetic example -----------------------------------------------
% Applying delay-cordinate DMD to a simple signal which is a linear 
%  span of sines and cosines with random frequencies between 0 and 20.  

% Generate data
dt = 0.1;
t = 0:dt:40;
numOsc = 2; 
x = zeros(1,length(t));

% Randomize frequencies and choice of sine or cosine
om = randi([0 20],numOsc,1);
sin_or_cos = randi([0 1],numOsc,1);
coeffs = 10*rand(numOsc,1) - 5;

% Create signal: x(t) = c_1*f_1(om_1*t) + c_2*f_2(om_2*t) + ... 
%   f_i are randomized to be either sine or cosine
for i = 1:numOsc
   
    if sin_or_cos(i) == 1
        x = x + coeffs(i)*sin(om(i)*t);
    else
        x = x + coeffs(i)*cos(om(i)*t);
    end
    
end

% Hankel matrix
delays = 2*numOsc; % tau = delays - 1
xd = hankel(x(1:delays+1),x(delays+1:end))';

% x_n+1 = a_1*x_{n - d + 1} + a_2*x_{n - d + 2} ... a_d*x_{n}  
xnp1 = xd(:,delays+1);
xn = xd(:,1:delays);
a = xn\xnp1;

% Create predicted signal
xpred = zeros(length(x),1);
xpred(1:delays) = x(1:delays); %xpred(1) = x(1), ... , xpred(d) = x(d) 

for j = delays+1:length(x)
    for k = 1:delays
       xpred(j) = xpred(j) + a(k)*xpred(j-delays+k-1); 
    end
end

% Plot results
plot(t,x,'b--','LineWidth',2)
hold on
plot(t,xpred,'r','LineWidth',1)
xlabel('t')
set(gca,'Fontsize',16,'Xlim',[0, 20])
legend('True Signal','Reconstructed Signal')

%% Generate Nonlinear Oscillation Example ---------------------------------
% Applying delay-cordinate DMD to the highly nonlinear Van der Pol signal.  

% Clean workspace
clear all; close all; clc

% Simulate ODE
dt = 0.05;
t = 0:dt:200;
x0 = [2; 2];
[t, x] = ode45(@(t,x) VdP(t,x),t,x0);

% Plot solution
subplot(2,1,1) % x(t)
plot(t,x(:,1),'b','Linewidth',2)
xlabel('t')
ylabel('x(t)')
set(gca,'Fontsize',16)
axis tight

subplot(2,1,2) % y(t)
plot(t,x(:,2),'r','Linewidth',2)
xlabel('t')
ylabel('y(t)')
set(gca,'Fontsize',16)
axis tight

%% Time-Delay Coordinate DMD on VdP ---------------------------------------

close all

% Hankel matrix
delays = 1000;
xd = hankel(x(1:delays+1,1),x(delays+1:end,1))';

% x_n+1 = a_1*x_{n - d + 1} + a_2*x_{n - d + 2} ... a_d*x_{n}  
xnp1 = xd(:,delays+1);
xn = xd(:,1:delays);
a = xn\xnp1;

% Create predicted signal
xpred = zeros(length(x(:,1)),1);
xpred(1:delays) = x(1:delays,1); %xpred(1) = x(1), ... , xpred(d) = x(d) 

for j = delays+1:length(x)
    for k = 1:delays
       xpred(j) = xpred(j) + a(k)*xpred(j-delays+k-1); 
    end
end

% Plot results
figure('Position', [100, 500, 700, 250])
plot(t,x(:,1),'b--','LineWidth',2)
hold on
plot(t,xpred,'r','LineWidth',2)
xlabel('$t$','interpreter','latex')
title(['$\tau =$ ', num2str(delays)],'interpreter','latex')
set(gca,'Fontsize',16,'Xlim',[0, t(end)/2],'Ylim',[-2.2, 2.2])
legend('True Signal','Reconstructed Signal')

%% SVD on Hankel Matrix ---------------------------------------------------

close all

% Hankel matrix
delays = 1000;
xd = hankel(x(1:delays,1),x(delays:end,1));

% Apply SVD
[U, S, V] = svd(xd); % SVD of delay matrix

% Plot SVD Results
figure('Position', [100, 500, 700, 500])
subplot(3,1,1) % Singular values
plot(diag(S)/max(diag(S)),'ko','Linewidth',2)
ylabel('$\sigma_j/\sigma_1$','interpreter','latex')
title(['$\tau =$ ', num2str(delays)],'interpreter','latex')
set(gca,'Fontsize',16,'Xlim',[0.9 min(delays+0.1,100+0.1)])

subplot(3,1,2) % Right-singular vectors
plot(t(1:end-delays+1),V(:,1),'r','Linewidth',2)
hold on
plot(t(1:end-delays+1),V(:,2),'b--','Linewidth',2)
xlabel('$t$','interpreter','latex')
set(gca,'Fontsize',16,'Xlim',[0 t(end-delays)])
legend('$v_1(t)$','$v_2(t)$','interpreter','latex')

subplot(3,1,3) % Right-singular vectors
plot(t(1:end-delays+1),V(:,3),'r','Linewidth',2)
hold on
plot(t(1:end-delays+1),V(:,4),'b--','Linewidth',2)
xlabel('$t$','interpreter','latex')
set(gca,'Fontsize',16,'Xlim',[0 t(end-delays)])
legend('$v_3(t)$','$v_4(t)$','interpreter','latex')

%% Time-Delay DMD on Low Rank Approximation of Hankel Matrix --------------

close all

% Find rank with 95% of energy
energy = 0;
energyTotal = sum(diag(S));
r = 0;
while energy <= 0.95
    r = r + 1;
    energy = energy + S(r,r)/energyTotal;
end

% Apply DMD to first r columns of V
X1 = V(1:end-1,1:r)';
X2 = V(2:end,1:r)'; 
[U2, S2, V2] = svd(X1,'econ');
A = U2'*X2*V2*diag(1./diag(S2));

% Plot eigenvalues of A
mu = eig(A);
plot(mu,'.','MarkerSize',10)
xlabel('Real Part','interpreter','latex')
ylabel('Imaginary Part','interpreter','latex')
title('Eigenvalues of the DMD Matrix','interpreter','latex')
set(gca,'Fontsize',16)


%% VdP Right-Hand-Side

function rhs = VdP(t,x)
    rhs = [x(2); -x(1) + 10*(1 - x(1)^2)*x(2)];
end
