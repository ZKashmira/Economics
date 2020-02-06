% Find initial level of consumption in the Neoclassical Growth Model
% Boundary problem - Need to find initial condition that satisfies Transversality
% Code generated for Ben Moll's Advanced Macroeconomics class at the LSE - Lecture 3

% Shooting Algorithm for Neoclassical Growth Model

% Clear environment
clear all; clc;

% Set parameter values
A = 1;
alpha = 0.3;
sigma = 2;
rho = 0.05;
delta = 0.05;
step = 0.1;
N = 700;

% Critical value to stop algorithm
crit = 10^(-6);

% Steady State capital
k_star = (alpha*A/(rho+delta))^(1/(1-alpha));
c_star = A*(k_star^alpha) - delta*k_star;

% Initial condition for capital
k0 = k_star/2;

% Get initial guess
c_upper_bound = A*(k0.^alpha) - delta.*k0;
c0 = c_upper_bound*rand(1);
k_old = k0;
c_old = c0;
c_lower_bound = 0;
k_new = k0;
diff = k_new - k_star;

for i=1:100
    
    % Evolution of c0
    for n=1:N
        % Get new value of k and c
        k_new = step .* (A .* (k_old.^alpha) - delta .* k_old - c_old) + k_old;
        c_new = step .* c_old .* 1./sigma .* (alpha .* A .* (k_old^(alpha-1)) - rho - delta) + c_old;
        if k_new <= 0| c_new <= 0
            break
        else
            k_old = k_new;
            c_old = c_new;
        end
    end
    
    diff = k_new - k_star;
    
    % Check if converged, if not, repeat
    if abs(diff) < crit
        disp('Algorithm converged')
        break
    elseif diff > 0
        c_lower_bound = c0;
        c0 = (c_upper_bound + c_lower_bound)./2;
        c_old = c0;
        k_old = k0;
    else
        c_upper_bound = c0;
        c0 = (c_upper_bound + c_lower_bound)./2;
        c_old = c0;
        k_old = k0;
    end
end

% Store all sequence
c_seq = 1:700;
k_seq = 1:700;
k_seq(1) = k0;
c_seq(1) = c0;
for i=2:N
    k_seq(i) = step .* (A .* (k_seq(i-1).^alpha) - delta .* k_seq(i-1) - c_seq(i-1)) + k_seq(i-1);
    c_seq(i) = step .* c_seq(i-1) .* 1./sigma .* (alpha .* A .* (k_seq(i-1)^(alpha-1)) - rho - delta) + c_seq(i-1);
end

% Plot
plot(k_seq, c_seq)

k = [0 : 0.01 : 8];
k_lom = A * (k.^alpha) - delta * k;
plot(k, k_lom, k_seq, c_seq, '--', k_seq(1), c_seq(1), '*')
xline(k_star, 'r')
title('Saddle path - Neoclassical Growth Model')
text(k_seq(1)-0.8,c_seq(1)-0.05,'Initial Condition','Color','red','FontSize',10)
xticks([0 k0 k_star])
xticklabels({'0','k0 = 2.402','k^* = 4.804'})
yticks([c0 c_star])
yticklabels({'c0 = 1.0045', 'c^* = 1.3611'})
