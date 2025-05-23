%% Clear workspace and close figures
clear all;
close all;
clc;

% Define the base path for your folders containing the data
dataPath = "C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Matlab\adaptivefilterMT\My_own_experiment\NHS data\speech";
folders = dir(fullfile(dataPath, "*"));

% Filter Parameters
M = 1; % Filter order (adjusted to match first code)

[primary, fs] = audioread("C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Matlab\adaptivefilterMT\4_Baseline_My_Own_Experiment\NHS data\speech\Primary.wav");  % Noisy signal with heartbeat   
[noise, ~] = audioread("C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Matlab\adaptivefilterMT\4_Baseline_My_Own_Experiment\NHS data\speech\Secondary.wav"); % Reference noise signal     
            
% % Calculate initial SNR for noisy signal (using reference noise)
% initial_SNR = 10 * log10(sum(noise.^2) / sum((primary - noise).^2));

% Pad noisy signal
N = length(primary); % Signal length
padded_signal = [zeros(M-1, 1); primary]; 

%% LMS Filter
mu_LMS = 0.6; % Step size for LMS (adjusted from first code)
w_LMS = zeros(M, 1); % Initialize filter weights
output_LMS = zeros(N, 1); % Filter output

for n = 1:N
    u_vect = padded_signal(n:n+M-1); % Current input vector
    % Calculate error signal (based on the noisy signal and filter output)
    e = noise(n) - w_LMS' * u_vect; % Error signal using the noisy signal as the desired response
    w_LMS = w_LMS + mu_LMS * e * u_vect; % Update weights
    output_LMS(n) = w_LMS' * u_vect; % Filtered output
end

filtered_signal_LMS = primary - output_LMS; 

%% NLMS Filter
mu_NLMS = 0.05; % Step size for NLMS (kept from first code)
w_NLMS = zeros(M, 1); % Initialize filter weights
output_NLMS = zeros(N, 1); % Filter output
Eps = 1e-5; % Stability constant (adjusted from first code)

for n = 1:N
    u_vect = padded_signal(n:n+M-1); % Current input vector
    mu_adapt = mu_NLMS / (Eps + norm(u_vect)^2); % Adaptive step size
    % Calculate error signal (based on the noisy signal and filter output)
    e = noise(n) - w_NLMS' * u_vect; % Error signal using the noisy signal as the desired response
    w_NLMS = w_NLMS + mu_adapt * e * u_vect; % Update weights
    output_NLMS(n) = w_NLMS' * u_vect; % Filtered output
end

filtered_signal_NLMS = primary - output_NLMS; 

%% RLS Filter
lambda = 0.95; % Forgetting factor  =  0.1666    =  1 - 1 / (0.1 * M)
delta = 0.01; % Initialization constant
P = (1 / delta) * eye(M); % Initialize inverse correlation matrix
w_RLS = zeros(M, 1); % Initialize filter weights
padded_signal = [sqrt(delta) * randn(M-1, 1); primary]; % Pad noisy signal
output_RLS = zeros(N, 1); % Filter output

for n = 1:N
    u_vect = padded_signal(n:n+M-1); % Current input vector
    PI = P * u_vect; % Intermediate calculation
    gain_k = PI / (lambda + u_vect' * PI); % Gain
    % Calculate error signal (based on the noisy signal and filter output)
    e = noise(n) - w_RLS' * u_vect; % Error signal using the noisy signal as the desired response
    w_RLS = w_RLS + gain_k * e; % Update weights
    P = P / lambda - gain_k * (u_vect' * P) / lambda; % Update P matrix
    output_RLS(n) = w_RLS' * u_vect; % Filtered output
end

filtered_signal_RLS = primary - output_RLS; 


% %% Calculate SNR Improvement (using the noisy signal as the clean signal)
% filtered_SNR_LMS = 10 * log10(sum(primary.^2) / sum((primary - filtered_signal_LMS).^2));
% filtered_SNR_NLMS = 10 * log10(sum(primary.^2) / sum((primary - filtered_signal_NLMS).^2));
% filtered_SNR_RLS = 10 * log10(sum(primary.^2) / sum((primary - filtered_signal_RLS).^2));
% 
% % Display Results
% fprintf('Initial SNR (before filtering): %.2f dB\n', initial_SNR);
% fprintf('SNR after LMS Filter: %.2f dB\n', filtered_SNR_LMS);
% fprintf('SNR after NLMS Filter: %.2f dB\n', filtered_SNR_NLMS);
% fprintf('SNR after RLS Filter: %.2f dB\n', filtered_SNR_RLS);

% Normalize filtered signals
filtered_signal_LMS = filtered_signal_LMS / max(abs(filtered_signal_LMS));
filtered_signal_NLMS = filtered_signal_NLMS / max(abs(filtered_signal_NLMS));
filtered_signal_RLS = filtered_signal_RLS / max(abs(filtered_signal_RLS));
primary = primary / max(abs(primary));  
noise = noise / max(abs(noise));  

%% Plot Mel Spectrogram for LMS, NLMS, and RLS in one figure

windowLength = round(0.14*fs);  % 140 ms window length
hop_size = round(0.02*fs);      % 20 ms hop size
overlap = windowLength - hop_size; 
numBands = 128;                 % Number of Mel bands
hannWin = hann(windowLength, 'periodic'); 

% Compute Mel spectrogram for primary signal
[s_u, f, t] = melSpectrogram(primary, fs, ...
                            "Window", kaiser(windowLength,18), ...
                            'OverlapLength', overlap, ...
                            'NumBands', numBands);

s_u_db = 10 * log10(s_u + eps);

% Compute Mel spectrogram for noise signal
[s_d, ~, ~] = melSpectrogram(noise, fs, ...
                            "Window", kaiser(windowLength,18), ...
                            'OverlapLength', overlap, ...
                            'NumBands', numBands);

s_d_db = 10 * log10(s_d + eps);

t = t - t(1); % Force time axis to start from 0

% Compute Mel spectrogram for filtered LMS signal
[s_filtered_LMS, ~, ~] = melSpectrogram(filtered_signal_LMS, fs, ...
                            "Window", kaiser(windowLength,18), ...
                            'OverlapLength', overlap, ...
                            'NumBands', numBands);

s_filtered_LMS_db = 10 * log10(s_filtered_LMS + eps);

% Compute Mel spectrogram for filtered NLMS signal
[s_filtered_NLMS, ~, ~] = melSpectrogram(filtered_signal_NLMS, fs, ...
                            "Window", kaiser(windowLength,18), ...
                            'OverlapLength', overlap, ...
                            'NumBands', numBands);

s_filtered_NLMS_db = 10 * log10(s_filtered_NLMS + eps);

% Compute Mel spectrogram for filtered RLS signal
[s_filtered_RLS, ~, ~] = melSpectrogram(filtered_signal_RLS, fs, ...
                            "Window", kaiser(windowLength,18), ...
                            'OverlapLength', overlap, ...
                            'NumBands', numBands);

s_filtered_RLS_db = 10 * log10(s_filtered_RLS + eps);

% Plot Mel spectrograms in a single figure (5 subplots)
figure;

% Primary signal
subplot(5,1,1); 
imagesc(t, f, s_u_db); 
axis xy; 
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Primary Signal');
colorbar;
colormap jet;
ylim([0 4000]);
xlim([0 5]);

% Noise signal
subplot(5,1,2);
imagesc(t, f, s_d_db); 
axis xy; 
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Noise Signal');
colorbar;
colormap jet;
ylim([0 4000]);
xlim([0 5]);

% Filtered LMS signal
subplot(5,1,3); 
imagesc(t, f, s_filtered_LMS_db); 
axis xy; 
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Filtered LMS Signal');
colorbar;
colormap jet;
ylim([0 4000]);
xlim([0 5]);

% Filtered NLMS signal
subplot(5,1,4);
imagesc(t, f, s_filtered_NLMS_db); 
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Filtered NLMS Signal');
colorbar;
colormap jet;
ylim([0 4000]);
xlim([0 5]);

% Filtered RLS signal
subplot(5,1,5);
imagesc(t, f, s_filtered_RLS_db); 
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Filtered RLS Signal');
colorbar;
colormap jet;
ylim([0 4000]);
xlim([0 5]);

% Set figure size and position
set(gcf, 'Units', 'inches', 'Position', [0, 0, 8.266666666666666, 9.866666666666667]);


% % Save the figure
% tightfig();
% saveas(gcf, fullfile(foldername, 'Mel_Spectrogram - My Heart.pdf'));


%% Visualization of Results
figure;
subplot(5, 1, 1);
plot(primary);
title('Heart + Noise signal');
xlabel('Sample Number');
ylabel('Amplitude');
xlim([0 length(primary)]);

subplot(5, 1, 2);
plot(noise);
title('Reference Noise Signal');
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

%% Save Filtered Signals
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
saveas(gcf, 'MOE_heartbeat_speech.pdf');

function tightfig()
    % Tighten the figure by removing excess whitespace
    set(gcf, 'Units', 'Inches');
    pos = get(gcf, 'Position');
    set(gcf, 'PaperUnits', 'Inches');
    set(gcf, 'PaperSize', [pos(3) pos(4)]);
    set(gcf, 'PaperPosition', [0 0 pos(3) pos(4)]);
end
