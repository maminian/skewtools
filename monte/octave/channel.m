% Main script for Monte Carlo for channel flow.
% Domain: [0,\infty) \times [0,1].
clear all; close all; clc

% Set parameters.
N = 10^2;      % number of walkers
T = 2;                  % final time
Pe = 2500;                  % Peclet number
dt = 10^-4;                % time step

t = 0:dt:T;

u = @(z) 4*Pe*(z.*(1-z)-1/6); % flow

% Define time and path vectors.


X = zeros(N,1);
Z = zeros(N,1);

tic

% Set initial conditions.

X = zeros(N,1);
%Z = 0.5*ones(N,1);
for i=0:N-1
     Z(i+1) = i/(N-1);
end

% Iterate time.
for k=2:length(t)

     % X straightforward enough because of simple Brownian motion.
     W_1 = sqrt(2)*randn(N,1);
     
     X = X + dt*u(Z) + 0*sqrt(dt)*W_1;
     
     % Z requires Brownian motion with reflective boundaries at 0,1.
     W_2 = sqrt(2)*randn(N,1);
     
     Z = Z + 0*sqrt(dt)*W_2;
     
     % Need to check if any of Z(:,k) lies outside [0,1]; if not,
     % we need to modify them relative to Wtemp to reflect the opposite
     % direction.
     for i=1:N
          if (Z(i) > 1)
               residual = Z(i) - 1;
               Z(i) = 1 - residual;
          elseif (Z(i) < 0)
               residual = Z(i) - 0;
               Z(i) = 0 - residual;
          end
     end
     
     channelvar(k) = var(X);
     channelskew(k) = skewness(X);
     
end


semilogx(t,channelskew)
xlabel('Time')
ylabel('Skewness')
mytitle = strcat('Skewness in channel with Pe=',num2str(Pe),', N=',num2str(N),', dt=',num2str(dt));
title(mytitle)

toc


