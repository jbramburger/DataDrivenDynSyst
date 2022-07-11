% -------------------------------------------------------------------------
% Neural Networks
%
% TBD
%
% Sigmoid activation function
%
% This script accompanies Section 3.2 of Computational Methods for
% Dynamical Systems. 
%
% Author: Jason J. Bramburger
% -------------------------------------------------------------------------

% Clean workspace
clear all; close all; clc

%% Initializations

% Network parameters
lr = 0.1;
n = 100; % hidden layer width

% Parameter initializations
w1 = rand(n,2);
w2 = rand(1,n);
b = zeros(n,1);

% Training data generation
    % To reproduce the example from the text load nn_data.mat
% load nn_data.mat
k = 1000;
x = 2*rand(2,k) - 1; % x is in R^2
y = dot(x,x); % y = x1^2 + x2^2

% Plot training data points
[x1p, x2p] = meshgrid(linspace(-1,1,1000),linspace(-1,1,1000));
figure(1)
hold on
view(-45,25)
surf(x1p,x2p,x1p.^2 + x2p.^2)
plot3(x(1,:),x(2,:),y,'.','Color',[1 69/255 79/255],'MarkerSize',25)
grid on
alpha 0.6
shading interp
colormap gray
xlabel('$x_1$','Interpreter','Latex')
ylabel('$x_2$','Interpreter','Latex')
zlabel('$y$','Interpreter','Latex')
set(gca,'FontSize',16)

%% Network training

% Number of gradient descent iterations
steps = 1e5;

% Initial loss
layer1 = ( 1./( 1 + exp(1).^(-w1*x - repmat(b,1,k)) ) );
network = w2*layer1;
loss(1) = (1/k)*norm(y - network)^2;

fprintf('step: 0 ----> loss %12.8f\n',0,loss(1));

% Gradient descent
for s = 1:steps
  
    % Compute gradients
     % Network gradients
    dgdw2 = layer1;
    dgdw1_1 = w2'.*layer1.*(1 - layer1).*x(1,:);
    dgdw1_2 = w2'.*layer1.*(1 - layer1).*x(2,:);
    dgdb = w2'.*layer1.*(1 - layer1);
    
     % Loss gradients
    dLdw2 = -(2/k)*dgdw2*(y - network)';
    dLdw1 = [-(2/k)*dgdw1_1*(y - network)' -(2/k)*dgdw1_2*(y - network)'];
    dLdb = -(2/k)*dgdb*(y - network)';
    
    % Updates
    b = b - lr*dLdb;
    w1 = w1 - lr*dLdw1;
    w2 = w2 - lr*dLdw2';
    
    % Compute new loss value
    layer1 = ( 1./( 1 + exp(1).^(-w1*x - repmat(b,1,k)) ) );
    network = w2*layer1;
    loss(s+1) = (1/k)*norm(y - network)^2;
    
    % Print loss
    fprintf('step: %12.0d ----> loss %12.8f\n',s+1,loss(s+1));
    
end

%% Plot Solution Surface

xspace = [x1p(:) x2p(:)]';
trainedNN =  w2*( 1./( 1 + exp(1).^(-w1*xspace - repmat(b,1,1000^2)) ) );
trainedNN = reshape(trainedNN,[1000 1000]);

figure(2)
hold on
view(-45,25)
surf(x1p,x2p,trainedNN)
%plot3(x(1,:),x(2,:),y,'.','Color',[1 69/255 79/255],'MarkerSize',25)
grid on
alpha 0.8
shading interp
colormap gray
xlabel('$x_1$','Interpreter','Latex')
ylabel('$x_2$','Interpreter','Latex')
zlabel('$y$','Interpreter','Latex')
set(gca,'FontSize',16)

%% Plot error

figure(3)
hold on
surf(x1p,x2p,abs(x1p.^2 + x2p.^2 - trainedNN))
plot3(x(1,:),x(2,:),y,'.','Color',[0.8 0.8 0.8],'MarkerSize',25)
view(0,90)
shading interp
colorbar
set(gca,'FontSize',16)




