function [x1UPO, x2UPO, x3UPO, period] = Gissinger_UPO(sequence)

% Function for producing UPOs in the flow of the chaotic Gissinger system. 
%
% Note: Sequences must be input as strings, meaning they are put in as '1'
% instad of 1, for example.
%
% Inputs:   - Sequences of 0's, 1's, and 2's 
%               - Acceptable sequences are listed in Table XX in Section
%               6.5 of textbook
%               - Default sequence is 1
%               - Invalid inputs are changed to the default
% 
% Outputs:  - x1UPO: the x1 component of the UPO
%           - x2UPO: the x2 component of the UPO
%           - x3UPO: the x3 component of the UPO
%           - period: the temporal period of the UPO
%
% This function accompanies Section 6.5 of Data-Driven Methods for
% Dynamic Systems. 
%
% Author: Jason J. Bramburger
% -------------------------------------------------------------------------

% Gissinger parameters
nu = 0.1;
gam = 0.85;
mu = 0.12;

% Integration parameters
m = 3; %Dimension of ODE
dt = 0.05;
tspan = 0:dt:20;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,m));

% Set default input
if nargin == 0 
    sequence = 1;
end

% Initial conditions in Poincare section
if strcmp(sequence,'1') == 1
    % 1 orbit periodic orbit
    x0(1,:) = [-1.2179209 1.2179209 1.3286662];
elseif strcmp(sequence,'01') == 1
    % 01 orbit periodic orbit
    x0(1,:) = [-1.073515 1.073515  1.1201103];
    x0(2,:) = [-1.2816105 1.2816105 1.4147233];
elseif strcmp(sequence,'12') == 1
    % 12 orbit periodic orbit
    x0(1,:) = [-1.1593777 1.1593777  1.247489];
    x0(2,:) = [-1.43768131 1.43768131 1.66258233];
elseif strcmp(sequence,'02') == 1
    % 02 orbit periodic orbit
    x0(1,:) = [-1.1006762 1.1006762 1.160484];
    x0(2,:) = [-1.4226933 1.4226933 1.6029335];
elseif strcmp(sequence,'001') == 1
    % 001 orbit periodic orbit
    x0(1,:) = [-1.0235752 1.0235752  1.039149];
    x0(2,:) = [-1.0818121 1.0818121 1.1322173];
    x0(3,:) = [-1.3238802 1.3238802 1.4743656];
elseif strcmp(sequence,'002') == 1
    % 002 orbit periodic orbit
    x0(1,:) = [-1.0275182 1.0275182  1.0456374];
    x0(2,:) = [-1.0889106 1.0889106 1.1433766];
    x0(3,:) = [-1.3570137 1.3570137 1.5177892];
elseif strcmp(sequence,'012') == 1
    % 012 orbit periodic orbit
    x0(1,:) = [-1.0456371 1.0456371 1.0748404];
    x0(2,:) = [-1.176175 1.176175  1.269483];
    x0(3,:) = [-1.3761501 1.3761501 1.543645];
elseif  strcmp(sequence,'011') == 1
    % 011 orbit periodic orbit
    x0(1,:) = [-1.0505573 1.0505573 1.0827899];
    x0(2,:) = [-1.1940813 1.1940813 1.2949812];
    x0(3,:) = [-1.2935671 1.2935671 1.4306426];
elseif  strcmp(sequence,'0012') == 1
    % 0012 orbit periodic orbit
    x0(1,:) = [-1.0217707 1.0217707 1.0364347];
    x0(2,:) = [-1.0475534 1.0475534 1.0779362];
    x0(3,:) = [-1.183928 1.183928  1.2805532];
    x0(4,:) = [-1.3434799 1.3434799 1.4995136];
elseif  strcmp(sequence,'0011') == 1
    % 0011 orbit periodic orbit
    x0(1,:) = [-1.0218817 1.0218817 1.036612];
    x0(2,:) = [-1.0484523 1.0484523 1.0793886];
    x0(3,:) = [-1.1874754 1.1874754 1.2856084];
    x0(4,:) = [-1.3288637 1.3288637  1.4792479];
elseif  strcmp(sequence,'0111') == 1
    % 0111 orbit periodic orbit
    x0(1,:) = [-1.058983 1.058983  1.0964093];
    x0(2,:) = [-1.2250694 1.2250694 1.3387289];
    x0(3,:) = [-1.1970693 1.1970693 1.2992151];
    x0(4,:) = [-1.2900715 1.2900715 1.4259255];
elseif  strcmp(sequence,'0112') == 1
    % 0112 orbit periodic orbit
    x0(1,:) = [-1.0614609 1.0614609 1.1004165];
    x0(2,:) = [-1.2328694 1.2328694 1.3496767];
    x0(3,:) = [-1.1733215 1.1733215 1.2654027];
    x0(4,:) = [-1.3878694 1.3878694 1.559321];
elseif  strcmp(sequence,'0102') == 1
    % 0102 orbit periodic orbit
    x0(1,:) = [-1.0704143 1.0704143 1.1150302];
    x0(2,:) = [-1.2683897 1.2683897 1.396983];
    x0(3,:) = [-1.0963032 1.0963032 1.1541272];
    x0(4,:) = [-1.3947937 1.3947937 1.5683061];
elseif  strcmp(sequence,'0212') == 1
    % 0212 orbit periodic orbit
    x0(1,:) = [-1.1100662 1.1100662 1.1741377];
    x0(2,:) = [-1.4736025 1.4736025 1.6684644];
    x0(3,:) = [-1.1651825 1.1651825 1.2537521];
    x0(4,:) = [-1.4313378 1.4313378 1.6140159];
elseif  strcmp(sequence,'0211') == 1
    % 0211 orbit periodic orbit
    x0(1,:) = [-1.1170244 1.1170244 1.1842479];
    x0(2,:) = [-1.498775 1.498775  1.7009403];
    x0(3,:) = [-1.2054536 1.2054536 1.3110778];
    x0(4,:) = [-1.2567233 1.2567233 1.3814785];
elseif  strcmp(sequence,'1211') == 1
    % 1211 orbit periodic orbit
    x0(1,:) = [-1.147692 1.147692  1.22858];
    x0(2,:) = [-1.5010884 1.5010884 1.7039268];
    x0(3,:) = [-1.2093925 1.2093925 1.3166412];
    x0(4,:) = [-1.2418877 1.2418877 1.3620665];
end

% Create initial guess
perMap = length(x0(:,1)); % period of UPO in Poincare map
xinit = [];
yinit = [];
zinit = [];

for p = 1:perMap

    % Initialize trajectory
    [~,sol] = ode45(@(t,x) Gissinger(x,mu,nu,gam),tspan,x0(p,:),options);

    ind = 0;
    for j = 1:length(sol(:,1))-1
       if  ((sol(j,1) + sol(j,2)) < 0) && ((sol(j+1,1) + sol(j+1,2)) >= 0)
            ind = j+1;

            xinit = [xinit; sol(1:ind,1)];
            yinit = [yinit; sol(1:ind,2)];
            zinit = [zinit; sol(1:ind,3)];

            break
        end 
    end

end

% Start clock
tic

% Root finding algorithm
T = dt*length(xinit(:,1)); % Initial guess for period
N = length(xinit(:,1));

% Time derivative
D = sparse(1:N-1,[2:N-1 N],ones(N-1,1)/2,N,N);
D(N,1) = 0.5; % Periodic BCs
D = (D - D')/dt;

% Add period to initial guess
uinit = [];
uinit = [xinit; yinit; zinit; T]; % Initial guess; % Initial guess

% fsolve options
options=optimset('Display','iter','Jacobian','on','MaxIter',10000,'Algorithm','levenberg-marquardt','TolFun',1e-15,'TolX',1e-15);

% call fsolve
xNew = fsolve(@(x) Periodic(x,D,uinit,mu,nu,gam,N),uinit,options);

toc %end timer

% Recover solutions (flip x & y because the conjugacy flips them in
% latent map)
x1UPO = xNew(1:N); 
x2UPO = xNew(N+1:2*N);
x3UPO = xNew(2*N+1:3*N);
period = T*xNew(end);

% Plot the solution in 3D
plot3([x1UPO; x1UPO],[x2UPO; x2UPO],[x3UPO; x3UPO],'Color',[36/255 122/255 254/255],'Linewidth',5)
set(gca,'fontsize',16)
xlabel('$x(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$y(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
zlabel('$z(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
grid on

end % end of Gissinger_UPO function

%% Gissinger right-hand-side

function dx = Gissinger(x,mu,nu,gamma)
    
    % Equilibria for scaling
    xstar = sqrt(nu + gamma*sqrt(nu/mu));
    ystar = sqrt(mu + gamma*sqrt(mu/nu));
    zstar = -sqrt(nu*mu) - gamma;
    
    % Rescaled variables
    x1hat = x(1)*xstar;
    x2hat = x(2)*ystar;
    x3hat = x(3)*zstar;

    dx = [(mu*x1hat - x2hat*(x3hat + gamma))/xstar; (-nu*x2hat + x1hat*(x3hat + gamma))/ystar; (-x3hat + x1hat*x2hat)/zstar];
    
end

%% Function whose roots are periodic orbits

function [F,J] = Periodic(xin,D,xinit,mu,nu,gamma,N)

  x = xin(1:N);
  y = xin(N+1:2*N);
  z = xin(2*N+1:3*N);
  T = xin(end);
  
  % Equilibria for scaling
  xstar = sqrt(nu + gamma*sqrt(nu/mu));
  ystar = sqrt(mu + gamma*sqrt(mu/nu));
  zstar = -sqrt(nu*mu) - gamma;

  % Right-hand side
  F(1:N) =  D*x - T*( (mu*x*xstar - ystar*y.*( z*zstar + gamma )) /xstar );
  F(N+1:2*N) =  D*y - T*( (-nu*y*ystar + x*xstar.*(z*zstar + gamma))/ystar );
  F(2*N+1:3*N) =  D*z - T*( (-z*zstar + xstar*ystar*x.*y)/zstar );
  F(3*N+1) =  (dot(D*xinit(1:N),xinit(1:N)-x) + dot(D*xinit(N+1:2*N),xinit(N+1:2*N)-y) + dot(D*xinit(2*N+1:3*N),xinit(2*N+1:3*N)-z));

  % Jacobian
  if nargout > 1
      e = ones(N,1);
      
      J = sparse(3*N+1,3*N+1);
      
      % First component
      J(1:N,1:N) = D + spdiags( -T*mu*e , 0, N, N);
      J(1:N,N+1:2*N) = spdiags( T*ystar*(z*zstar + gamma)/xstar, 0, N, N);
      J(1:N,2*N+1:3*N) = spdiags( T*(ystar*zstar/xstar)*y , 0, N, N);
      J(1:N,3*N+1) = -(mu*x*xstar - ystar*y.*(z*zstar + gamma))/xstar; 
      
      % Second component 
      J(N+1:2*N,1:N) = spdiags( -T*xstar*(z*zstar + gamma)/ystar , 0, N, N);
      J(N+1:2*N,N+1:2*N) = D + spdiags( T*nu*e , 0, N, N);
      J(N+1:2*N,2*N+1:3*N) = spdiags( -T*x*xstar*zstar/ystar , 0, N, N); 
      J(N+1:2*N,3*N+1) = -(-nu*y*ystar + x*xstar.*(z*zstar + gamma))/ystar; 
      
      % Third component
      J(2*N+1:3*N,1:N) = spdiags( -T*xstar*ystar*y/zstar , 0, N, N);
      J(2*N+1:3*N,N+1:2*N) = spdiags( -T*xstar*ystar*x/zstar , 0, N, N);
      J(2*N+1:3*N,2*N+1:3*N) = D + spdiags( T*e , 0, N, N);
      J(2*N+1:3*N,3*N+1) = -( (-z*zstar + xstar*ystar*x.*y)/zstar );
      
      % Fourth component
      J(3*N+1,1:N) = -D*xinit(1:N);
      J(3*N+1,N+1:2*N) = -D*xinit(N+1:2*N);
      J(3*N+1,2*N+1:3*N) = -D*xinit(2*N+1:3*N);
  end

end
