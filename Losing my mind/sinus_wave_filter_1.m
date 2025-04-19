clear all;
close all;
clc;

%% Simulate Signal and Noise
fs = 8000; % Sampling frequency
duration = 15; % Duration of the signal (in seconds)
t = 0:1/fs:duration; % Time vector
clean = sin(2 * pi * 10 * t)'; % Simulated clean signal (10 Hz sine wave)
noise = 0.5 * randn(size(t))'; % Simulated Gaussian noise
primary = clean + noise; % Combine signal and noise

% Calculate initial SNR
initial_SNR = snr(clean,noise);

% Filter Parameters
N = length(primary); % Signal length
M = 12; % Filter taps/order (adjusted to match first code)

%% LMS Filter
mu_LMS = 0.0840; % Step size for LMS (adjusted from first code) %% 0.0006
w_LMS = zeros(M, 1); % Initialize filter weights
padded_signal = [zeros(M-1, 1); primary]; % Pad noisy signal
output_LMS = zeros(N, 1); % Filter output

for n = 1:N
    u_vect = padded_signal(n:n+M-1); % Current input vector
    e = noise(n) - w_LMS' * u_vect; % Error signal
    w_LMS = w_LMS + mu_LMS * e * u_vect; % Update weights
    output_LMS(n) = w_LMS' * u_vect; % Filtered output
end 

filtered_signal_LMS = primary - output_LMS; 

%% NLMS Filter
mu_NLMS = 1; % Step size for NLMS (kept from first code)
w_NLMS = zeros(M, 1); % Initialize filter weights
padded_signal = [zeros(M-1, 1); primary]; % Pad noisy signal
output_NLMS = zeros(N, 1); % Filter output
Eps = 0.0001; % Stability constant (adjusted from first code)

for n = 1:N
    u_vect = padded_signal(n:n+M-1); % Current input vector
    mu_adapt = mu_NLMS / (Eps + norm(u_vect)^2); % Adaptive step size
    e = noise(n) - w_NLMS' * u_vect; % Error signal
    w_NLMS = w_NLMS + mu_adapt * e * u_vect; % Update weights
    output_NLMS(n) = w_NLMS' * u_vect; % Filtered output
end

filtered_signal_NLMS = primary - output_NLMS; 

%% RLS Filter
lambda = 1 - 1 / (0.1 * M); % Forgetting factor (adjusted from first code)
delta = 0.01; % Initialization constant (adjusted from first code)
P = (1 / delta) * eye(M); % Initialize inverse correlation matrix
w_RLS = zeros(M, 1); % Initialize filter weights
padded_signal = [sqrt(delta) * randn(M-1, 1); primary]; % Pad noisy signal
output_RLS = zeros(N, 1); % Filter output

for n = 1:N
    u_vect = padded_signal(n:n+M-1); % Current input vector
    PI = P * u_vect; % Intermediate calculation
    gain_k = PI / (lambda + u_vect' * PI); % Gain
    e = noise(n) - w_RLS' * u_vect; % Error signal
    w_RLS = w_RLS + gain_k * e; % Update weights
    P = P / lambda - gain_k * (u_vect' * P) / lambda; % Update P matrix
    output_RLS(n) = w_RLS' * u_vect; % Filtered output
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

%% Export filtered signals

% Normalize filtered signals    (To prevent "Warning: Data clipped when writing file") 
output_LMS_normalized = filtered_signal_LMS / max(abs(filtered_signal_LMS));
output_NLMS_normalized = filtered_signal_NLMS / max(abs(filtered_signal_NLMS));
output_RLS_normalized = filtered_signal_RLS / max(abs(filtered_signal_RLS));

% Get the current directory where the MATLAB script is located
currentFolder = fileparts(mfilename('fullpath'));

% Define the name of the new folder
newFolder = fullfile(currentFolder, 'Filtered_Audio_Files');

% Create the folder if it does not exist
if ~exist(newFolder, 'dir')
    mkdir(newFolder);
end

% Save Filtered Signals in the newly created folder
audiowrite(fullfile(newFolder, 'Filtered_LMS_1.wav'), output_LMS_normalized, fs);
audiowrite(fullfile(newFolder, 'Filtered_NLMS_1.wav'), output_NLMS_normalized, fs);
audiowrite(fullfile(newFolder, 'Filtered_RLS_1.wav'), output_RLS_normalized, fs);

fprintf('Filtered signals saved as audio files.\n\n');

% Tuning Parementers for tables
fprintf('Sampling rate: %f Hz\n', fs);
fprintf('Filter order (LMS, NLMS, RLS): %f \n', M);
fprintf('Step size (LMS): %f \n', mu_LMS);
fprintf('Step size (NLMS): %f \n', mu_NLMS);
fprintf('Forgetting factor (RLS): %f \n', lambda);


% Crop the figure and save as PDF
tightfig();
saveas(gcf, 'sinus_wave_filter_1.pdf');

function tightfig()
    % Tighten the figure by removing excess whitespace
    set(gcf, 'Units', 'Inches');
    pos = get(gcf, 'Position');
    set(gcf, 'PaperUnits', 'Inches');
    set(gcf, 'PaperSize', [pos(3) pos(4)]);
    set(gcf, 'PaperPosition', [0 0 pos(3) pos(4)]);
end
