%% Clear workspace and close figures
clear all;
close all;
clc;

%% Define the path for the data
%[primary, fs] = audioread("C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Matlab\adaptivefilterMT\3_Baseline_Company_Data\Data_ANC\Experiment_Data\Artifacts\NHS\1\primary.wav");   %noise + clean signal
%[noise, ~] = audioread("C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Matlab\adaptivefilterMT\3_Baseline_Company_Data\Data_ANC\Experiment_Data\Artifacts\NHS\1\secondary.wav");
%[clean, ~] = audioread("C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Matlab\adaptivefilterMT\3_Baseline_Company_Data\Data_ANC\Experiment_Data\Artifacts\NHS\1\ZCH0048.wav");

[primary, fs] = audioread("C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Data from AI heathway\Data_ANC\Experiment_Data\Hospital Ambient Noises\NHS\1\primary.wav");   %noise + clean signal
[noise, ~] = audioread("C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Data from AI heathway\Data_ANC\Experiment_Data\Hospital Ambient Noises\NHS\1\secondary.wav");
[clean, ~] = audioread("C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Data from AI heathway\Data_ANC\Experiment_Data\Hospital Ambient Noises\NHS\1\ZCH0019.wav");

% Find the minimum length
minLength = min([length(primary), length(noise), length(clean)]);

% Truncate signals to the same length
primary = primary(1:minLength);
noise = noise(1:minLength);
clean = clean(1:minLength);

% Calculate initial SNR for noisy signal
initial_SNR = 10 * log10(sum(clean.^2) / sum((clean-noise).^2));
fprintf('Initial SNR (before filtering): %.2f dB\n', initial_SNR);

% Filter Parameters
M = 12; % Filter order (adjusted to match first code)
N = length(primary); % Signal length

%% LMS Filter
mu_LMS = 0.1; % Step size for LMS (adjusted from first code)
w_LMS = zeros(M, 1); % Initialize filter weights
padded_signal = [zeros(M-1, 1); primary]; % Pad input signal
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
padded_signal = [zeros(M-1, 1); primary]; % Pad input signal
output_NLMS = zeros(N, 1); % % Initialize output
Eps = 0.0001; % Small constant for stability

for n = 1:N
    u_vect = padded_signal(n:n+M-1); % Current input vector
    mu_adapt = mu_NLMS / (Eps + norm(u_vect)^2); % Adaptive step size
    e = noise(n) - w_NLMS' * u_vect; % Error signal
    w_NLMS = w_NLMS + mu_adapt * e * u_vect; % Update weights
    output_NLMS(n) = w_NLMS' * u_vect; % Filtered output
end

filtered_signal_NLMS = primary - output_NLMS; 

%% RLS Filter
lambda = 1 - 1 / (0.1 * M); % Forgetting factor  =  0.1666 
delta = 0.01; % Small constant for initialization
P = 1 / delta * eye(M); % Initialize inverse correlation matrix
w_RLS = zeros(M, 1); % Initialize filter weights
padded_signal = [sqrt(delta) * randn(M-1, 1); primary]; % Pad input signal
output_RLS = zeros(N, 1); % Initialize output

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

% filtered_SNR_LMS = snr(clean, filtered_signal_LMS - clean);
% filtered_SNR_NLMS = snr(clean, filtered_signal_NLMS - clean);
% filtered_SNR_RLS = snr(clean, filtered_signal_RLS - clean);

% Display Results
fprintf('SNR after LMS Filter: %.2f dB\n', mean(filtered_SNR_LMS));
fprintf('SNR after NLMS Filter: %.2f dB\n', mean(filtered_SNR_NLMS));
fprintf('SNR after RLS Filter: %.2f dB\n', mean(filtered_SNR_RLS));

%% Visualization of Results
figure;
subplot(5, 1, 1);
plot(clean);
title('Clean Signal');
xlabel('Sample Number');
ylabel('Amplitude');
xlim([0 length(primary)]);

subplot(5, 1, 2);
plot(primary);
title('Noisy Signal');
xlabel('Sample Number');
ylabel('Amplitude');
xlim([0 length(primary)]);

subplot(5, 1, 3);
plot(filtered_signal_LMS);
title('Filtered Signal (LMS)');
xlabel('Sample Number');
ylabel('Amplitude');
xlim([0 length(primary)]);

subplot(5, 1, 4);
plot(filtered_signal_NLMS);
title('Filtered Signal (NLMS)');
xlabel('Sample Number');
ylabel('Amplitude');
xlim([0 length(primary)]);

subplot(5, 1, 5);
plot(filtered_signal_RLS);
title('Filtered Signal (RLS)');
xlabel('Sample Number');
ylabel('Amplitude');
xlim([0 length(primary)]);

% figure
% plot(primary)
% hold on
% plot(filtered_signal_LMS,'r--')
% hold off
% title('Primary vs filtered LMS');
% legend('Received','Filtered')
% grid on
% ylabel('Amplitude')
% xlabel('Time (n)')

%% Export filtered signals

% Normalize filtered signals
filtered_signal_LMS = filtered_signal_LMS / max(abs(filtered_signal_LMS));
filtered_signal_NLMS = filtered_signal_NLMS / max(abs(filtered_signal_NLMS));
filtered_signal_RLS = filtered_signal_RLS / max(abs(filtered_signal_RLS));

% Save Filtered Signals
audiowrite('Filtered_LMS.wav', filtered_signal_LMS, fs);
audiowrite('Filtered_NLMS.wav', filtered_signal_NLMS, fs);
audiowrite('Filtered_RLS.wav', filtered_signal_RLS, fs);

fprintf('Filtered signals saved as audio files.\n\n');

% Tuning Parementers for tables
fprintf('Sampling rate: %f Hz\n', fs);
fprintf('Filter order (LMS, NLMS, RLS): %f \n', M);
fprintf('Step size (LMS): %f \n', mu_LMS);
fprintf('Step size (NLMS): %f \n', mu_NLMS);
fprintf('Forgetting factor (RLS): %f \n', lambda);

% Crop the figure and save as PDF
tightfig();
saveas(gcf, 'Hospital Ambient Noises-NHS-1.pdf');

function tightfig()
    % Tighten the figure by removing excess whitespace
    set(gcf, 'Units', 'Inches');
    pos = get(gcf, 'Position');
    set(gcf, 'PaperUnits', 'Inches');
    set(gcf, 'PaperSize', [pos(3) pos(4)]);
    set(gcf, 'PaperPosition', [0 0 pos(3) pos(4)]);
end
