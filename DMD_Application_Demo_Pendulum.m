clear all; 
close all;
clc;
%% Define time and space discretizations
%divides the space up into 180 points including endpoints (179 intervals)
t = linspace(0,0.5,10); 
% dt = time step advancing X1 to X2 (X to X')
dt = t(2) - t(1); 

%% Combine signals and make data matrix
theta_ic = [0.5; 0];       %initial conditions: theta(t=0)=0.5; dtheta(t=0)=0.
tspan = [0,10];
[t, theta] = ode45(@odeFun, tspan, theta_ic);

%% Create DMD data matrices
X=theta';
X1 = X(:, 1:end-1);        %X
X2 = X(:, 2:end);          %X'

%% SVD and rank-2 truncation
r = 2; % rank truncation r = min{m, n},
[U, S, V] = svd(X1, 'econ');
Ur = U(:, 1:r); 
Sr = S(1:r, 1:r);          %S^-1=Sigma inverse
Vr = V(:, 1:r);

%% Build Atilde and DMD Modes
Atilde = Ur'*X2*Vr/Sr;
[W, D] = eig(Atilde);      %generate e-values of Atilda using eig command
Phi = X2*Vr/Sr*W;          % DMD Modes

%% DMD Spectra
lambda = diag(D);          %generate e-vectors of Atilda using lambda command
omega = log(lambda)/dt;    %mapping e-values

%% Compute DMD Solution
x1 = X(:, 1);
b = Phi\x1;
time_dynamics = zeros(r,length(t)-1);
for iter = 1:length(t)
    time_dynamics(:,iter) = b.*exp(omega*t(iter));
end

X_dmd = Phi*time_dynamics;

%% Plot Functions
figure;
subplot(3,2,1);
plot(t,theta(:,1),'linewidth', 2); 
legend('X1', 'FontSize', 12);

subplot(3,2,2);
plot(t,theta(:,2),'linewidth', 2);
legend('X2', 'FontSize', 12);


subplot(3,2,3);
plot(t, theta,'linewidth', 2);
legend('X', 'FontSize', 12);

subplot(3,2,4); 
plot(t,real(X_dmd)','linewidth', 2);
legend('DMD-X1', 'DMD-X2','FontSize', 12);

subplot(3,2,5);
plot(t,real(X_dmd(1,:))',t,theta(:,1),'linewidth', 2);legend('theta', 'FontSize', 12);
legend('X1','DMD-X1', 'FontSize', 12);

subplot(3,2,6); 
plot(t,real(X_dmd(2,:))',t,theta(:,2),'linewidth', 2);
legend('X2','DMD-X2', 'FontSize', 12);


function dtheta = odeFun(~, theta)
   g = 9.8;
   l = 1;
   dtheta = zeros(2, 1);
   dtheta(1) = theta(2);          % =theta'
   dtheta(2) = -g/l*theta(1);     % =-g/l*theta
end
