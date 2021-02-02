%% Define time and space discretizations
xi = linspace(0,10,40) ;
t = linspace(0,pi,20);
dt = t(2) - t(1);
[Xgrid ,T] = meshgrid(xi,t);
%% Combine signals and make data matrix
theta_ic = [0.5; 0]; % initial conditions: theta(t=0)=0.5; dtheta(t=0)=0.
tspan = [0,10];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t, thteta] = ode45(@odeFun, tspan, theta_ic);
%% Create DMD data matrices
X=theta;
X1 = X(:, 1:end-1);
X2 = X(:, 2:end);

%% SVD and rank-1 truncation
r = 1; % rank truncation
[U, S, V] = svd(X1, 'econ');
Ur = U(:, 1:r); 
Sr = S(1:r, 1:r);
Vr = V(:, 1:r);

%% Build Atilde and DMD Modes
Atilde = Ur'*X2*Vr/Sr;
[W, D] = eig(Atilde);
Phi = X2*Vr/Sr*W;          % DMD Modes

%% DMD Spectra
lambda = diag(D);
omega = log(lambda)/dt;

%% Compute DMD Solution
x1 = X(:, 1);
b = Phi\x1;
time_dynamics = zeros(r,length(t));
for iter = 1:length(t)
    time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
end

X_dmd = Phi*time_dynamics;

figure;
subplot(2,2,1);
plot(t,theta(:,1)); axis equal

subplot(2,2,2);
plot(t,theta(:,2),'r'); axis equal

subplot(2,2,3);
plot(t, theta); axis equal


subplot(2,2,4); 
plot(real(X_dmd')); 

function dtheta = odeFun(~, theta)

   g = 9.8;
   l = 1;
   % theta(1) = theta, theta(2) = dtheta
  
   dtheta = zeros(2, 1);
   dtheta(1) = theta(2);         % =theta'
   dtheta(2) = -g/l*theta(1);    % =-g/l*theta
   
end
