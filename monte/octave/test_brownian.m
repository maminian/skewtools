% To demonstrate that sampling
% Brownian motion requires,
% for timestep dt, iid normal
% steps of variance dt, effected
% by multiplying by sqrt(dt).
%
% Any other scaling will give off
% a visual cue; the order of the
% values at t=5 will not stay
% in the same order for any scaling
% other than sqrt(dt).

clear all; close all; clc
n = 5; % Test time steps.

for k=1:n
     figure; hold on
     t = linspace(0,5,10^k);
     dt = t(2) - t(1);
               
     for derp=1:10
          X = randn(size(t));
          B = cumsum(sqrt(dt)*X);

          plot(t,B)
     end
     hold off
end
