clear all;
close all;
clc;

%% Simulate Signal and Noise
fs = 8000; % Sampling frequency
duration = 10; % Duration of the signal (in seconds)
t = 0:1/fs:duration; % Time vector
clean = sin(2 * pi * 10 * t)'; % Simulated clean signal (10 Hz sine wave)
noise = 0.5 * randn(size(t))'; % Simulated Gaussian noise
primary = clean + noise; % Combine signal and noise

% Calculate initial SNR
initial_SNR = snr(clean,noise);

% Filter Parameters
N = length(primary); % Signal length
M = 12; % Filter taps/order (adjusted to match first code)
mu_LMS = 0.0840; % Step size for LMS
mu_NLMS = 1; % Step size for NLMS
lambda = 1 - 1 / (0.1 * M); % Forgetting factor (adjusted from first code)


%% LMS Filter
w_LMS = zeros(M, 1); % Initialize filter weights
padded_signal = [zeros(M-1, 1); primary]; % Pad noisy signal
output_LMS = zeros(N, 1); % Filter output
mse_LMS = zeros(N, 1);
e_LMS = zeros(N,1);

for n = 1:N
    u_vect = padded_signal(n:n+M-1); % Current input vector
    output_LMS(n) = w_LMS' * u_vect; % Filtered output
    e_LMS(n) = noise(n) - output_LMS(n) ; % Error signal
    w_LMS = w_LMS + mu_LMS * e_LMS(n) * u_vect; % Update weights
    
    mse_LMS(n) = norm(mse_LMS(n)+(abs(e_LMS(n))^2));
end 

filtered_signal_LMS = primary - output_LMS; 

%% NLMS Filter
w_NLMS = zeros(M, 1); % Initialize filter weights
padded_signal = [zeros(M-1, 1); primary]; % Pad noisy signal
output_NLMS = zeros(N, 1); % Filter output
Eps = 0.0001; % Stability constant (adjusted from first code)
mse_NLMS = zeros(N, 1);
e_NLMS = zeros(N,1);


for n = 1:N
    u_vect = padded_signal(n:n+M-1); % Current input vector
    output_NLMS(n) = w_NLMS' * u_vect; % Filtered output
    e_NLMS(n) = noise(n) - output_NLMS(n); % Error signal
    mu_adapt = mu_NLMS / (Eps + norm(u_vect)^2); % Adaptive step size
   
    w_NLMS = w_NLMS + mu_adapt * e_NLMS(n) * u_vect; % Update weights
    mse_NLMS(n) = norm(mse_NLMS(n)+(abs(e_NLMS(n))^2));
end

filtered_signal_NLMS = primary - output_NLMS; 

%% RLS Filter
delta = 0.01; % Initialization constant (adjusted from first code)
P = (1 / delta) * eye(M); % Initialize inverse correlation matrix
w_RLS = zeros(M, 1); % Initialize filter weights
padded_signal = [sqrt(delta) * randn(M-1, 1); primary]; % Pad noisy signal
output_RLS = zeros(N, 1); % Filter output
mse_RLS = zeros(N, 1);
e_RLS = zeros(N,1);

for n = 1:N
    u_vect = padded_signal(n:n+M-1); % Current input vector
    output_RLS(n) = w_RLS' * u_vect; % Filtered output
    e_RLS(n) = noise(n) - output_RLS(n); % Error signal
    PI = P * u_vect; % Intermediate calculation
    gain_k = PI / (lambda + u_vect' * PI); % Gain
    w_RLS = w_RLS + gain_k * e_RLS(n); % Update weights
    P = P / lambda - gain_k * (u_vect' * P) / lambda; % Update P matrix
    
    mse_RLS(n) = norm(mse_RLS(n)+(abs(e_RLS(n))^2));
end

filtered_signal_RLS = primary - output_RLS; 

%% Calculate SNR Improvement
filtered_SNR_LMS = 10 * log10(sum(clean.^2) / sum((clean - filtered_signal_LMS).^2));
filtered_SNR_NLMS = 10 * log10(sum(clean.^2) / sum((clean - filtered_signal_NLMS).^2));
filtered_SNR_RLS = 10 * log10(sum(clean.^2) / sum((clean - filtered_signal_RLS).^2));

fprintf('Initial SNR: %.2f dB\n', initial_SNR);
fprintf('SNR after LMS Filter: %.2f dB\n', filtered_SNR_LMS);
fprintf('SNR after NLMS Filter: %.2f dB\n', filtered_SNR_NLMS);
fprintf('SNR after RLS Filter: %.2f dB\n', filtered_SNR_RLS);

%% Visualization of Results
figure;
subplot(5, 1, 1);
plot(clean);
title('Clean Signal');
xlabel('Sample Number');
ylabel('Amplitude');
ylim([-1.1 1.1]);
xlim([0 15000]);

subplot(5, 1, 2);
plot(primary);
title('Noisy Signal');
xlabel('Sample Number');
ylabel('Amplitude');
ylim([-1.1 1.1]);
xlim([0 15000]);

subplot(5, 1, 3); 
plot(clean-filtered_signal_LMS);
title('Error signal (LMS)');
xlabel('Sample Number');
ylabel('Amplitude');
ylim([-1.1 1.1]);
xlim([0 15000]);

subplot(5, 1, 4);
plot(clean-filtered_signal_NLMS);
title('Error signal (NLMS)');
xlabel('Sample Number');
ylabel('Amplitude');
ylim([-1.1 1.1]);
xlim([0 15000]);

subplot(5, 1, 5);
plot(clean-filtered_signal_RLS);
title('Error signal (RLS)');
xlabel('Sample Number');
ylabel('Amplitude');
ylim([-1.1 1.1]);
xlim([0 15000]);

%% Plot learning curves
figure(10)
plot(mse_LMS(1:1000),':k')
ylabel('ensemble-average squared error')
xlabel('number of iterations')
title('LMS - Convergence rate ')
axis([0 1000 0 1])
grid on
figure(11)
plot(mse_NLMS(1:1000),':k')
ylabel('ensemble-average squared error')
xlabel('number of iterations')
title('NLMS- Convergence rate ')
axis([0 1000 0 1])
grid on
figure(12)
plot(mse_RLS(1:1000),':k')
ylabel('ensemble-average squared error')
xlabel('number of iterations')
title('RLS - Convergence rate ')
axis([0 1000 0 1])
grid on
%% Export filtered signals

% Normalize filtered signals    (To prevent "Warning: Data clipped when writing file") 
output_LMS_normalized = filtered_signal_LMS / max(abs(filtered_signal_LMS));
output_NLMS_normalized = filtered_signal_NLMS / max(abs(filtered_signal_NLMS));
output_RLS_normalized = filtered_signal_RLS / max(abs(filtered_signal_RLS));

% Save Filtered Signals
audiowrite('Filtered_LMS.wav', output_LMS_normalized, fs);
audiowrite('Filtered_NLMS.wav', output_NLMS_normalized, fs);
audiowrite('Filtered_RLS.wav', output_RLS_normalized, fs);

fprintf('Filtered signals saved as audio files.\n\n');

% Tuning Parementers for tables
fprintf('Sampling rate: %f Hz\n', fs);
fprintf('Filter order (LMS, NLMS, RLS): %f \n', M);
fprintf('Step size (LMS): %f \n', mu_LMS);
fprintf('Step size (NLMS): %f \n', mu_NLMS);
fprintf('Forgetting factor (RLS): %f \n', lambda);


% Crop the figure and save as PDF
tightfig();
saveas(gcf, 'baseline_sim.pdf');

function tightfig()
    % Tighten the figure by removing excess whitespace
    set(gcf, 'Units', 'Inches');
    pos = get(gcf, 'Position');
    set(gcf, 'PaperUnits', 'Inches');
    set(gcf, 'PaperSize', [pos(3) pos(4)]);
    set(gcf, 'PaperPosition', [0 0 pos(3) pos(4)]);
end
