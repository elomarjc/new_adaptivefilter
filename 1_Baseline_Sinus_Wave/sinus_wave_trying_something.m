clear all;
close all;

%% Simulate Signal and Noise
fs = 8000; % Sampling frequency
duration = 15; % seconds
t = 0:1/fs:duration-1/fs; % time vector
clean = sin(2 * pi * 10 * t)'; % Clean signal
noise = 0.5 * randn(size(t))'; % Noise
primary = clean + noise;       % Observed signal

%% Adaptive filter parameters
N = length(t); % number of samples
R = 100;       % number of ensemble runs
nord = 25;     % filter order

% LMS
mu = 0.0500;

% NLMS
beta = 0.5400;

% RLS
delta = 0.001;
lambda = 1;

% Error performance initialization
MSE_LMS = zeros(N,1);
MSE_NLMS = zeros(N,1);
MSE_RLS = zeros(N,1);

%% Ensemble loop
for r = 1:R
    d = noise;           % Desired signal
    x = primary;         % Noisy signal
    v1 = noise;          % Noise component

    % Delayed reference signal
    n0 = 25; % delay
    x_del = zeros(N,1);
    x_del(1:end-n0) = x(n0+1:end);

    % Filter weights
    W_LMS = zeros(nord,1);
    W_NLMS = zeros(nord,1);
    W_RLS = zeros(nord,1);

    % RLS matrix
    P = (1/delta) * eye(nord);

    % Input vector
    U = zeros(nord,1);

    for i = 1:N
        % Update input buffer
        U = [x_del(i); U(1:end-1)];
        x_n = x(i);

        %% LMS
        y_LMS = W_LMS' * U;
        e_LMS = x_n - y_LMS;
        W_LMS = W_LMS + mu * e_LMS * U;
        MSE_LMS(i) = MSE_LMS(i) + abs(e_LMS)^2;

        %% NLMS
        y_NLMS = W_NLMS' * U;
        e_NLMS = x_n - y_NLMS;
        normU = norm(U)^2 + 1e-6; % Avoid divide-by-zero
        W_NLMS = W_NLMS + (beta / normU) * e_NLMS * U;
        MSE_NLMS(i) = MSE_NLMS(i) + abs(e_NLMS)^2;

        %% RLS
        g = (1/lambda) * P * U / (1 + (1/lambda) * U' * P * U);
        y_RLS = W_RLS' * U;
        e_RLS = x_n - y_RLS;
        W_RLS = W_RLS + g * e_RLS;
        P = (1/lambda) * (P - g * U' * P);
        MSE_RLS(i) = MSE_RLS(i) + abs(e_RLS)^2;
    end
end

%% Final MSE averaging
MSE_LMS = MSE_LMS / R;
MSE_NLMS = MSE_NLMS / R;
MSE_RLS = MSE_RLS / R;

%% Plot learning curves
figure;
plot(MSE_LMS(1:1000), 'r')
ylabel('MSE')
xlabel('Iterations')
title('LMS - Convergence')
axis([0 1000 0 1])
grid on

figure;
plot(MSE_NLMS(1:1000), 'b')
ylabel('MSE')
xlabel('Iterations')
title('NLMS - Convergence')
axis([0 1000 0 1])
grid on

figure;
plot(MSE_RLS(1:1000), 'g')
ylabel('MSE')
xlabel('Iterations')
title('RLS - Convergence')
axis([0 1000 0 1])
grid on
