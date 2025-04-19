%% Clear workspace and close figures
clear all;
close all;
clc;

%% Define the path for the data
[primary, fs] = audioread("C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Data from AI heathway\Data_ANC\Experiment_Data\Hospital Ambient Noises\NHS\1\primary.wav");   %noise + clean signal
[noise, ~] = audioread("C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Data from AI heathway\Data_ANC\Experiment_Data\Hospital Ambient Noises\NHS\1\secondary.wav");
[clean, ~] = audioread("C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Data from AI heathway\Data_ANC\Experiment_Data\Hospital Ambient Noises\NHS\1\ZCH0019.wav");

%% Define start time for trimming (1.5 seconds)
start_time = 1.7;  % in seconds
start_sample = round(start_time * fs);  % Convert time to sample index

%% Trim the signals to start from 1.5 seconds
primary = primary(start_sample:end);
noise = noise(start_sample:end);
clean = clean(start_sample:end);

% Find the minimum length
minLength = min([length(primary), length(noise), length(clean)]);

% Truncate signals to the same length
primary = primary(1:minLength);
noise = noise(1:minLength);
clean = clean(1:minLength);

%% Calculate initial SNR for audio files
initial_SNR = 10 * log10(sum(clean.^2) / sum((clean - noise).^2));

fprintf('Initial SNR (before any processing): %.2f dB\n', initial_SNR);
%% Fix #1: Bandpass Filtering to Remove Distortions
lowCutoff = 5; % Low cut frequency in Hz
highCutoff = 800; % High cut frequency in Hz
[b, a] = butter(4, [lowCutoff, highCutoff] / (fs/2), 'bandpass');
primary = filtfilt(b, a, primary);
noise = filtfilt(b, a, noise);

% SNR after fix #1 (Bandpass Filtering)
SNR_after_bandpass = 10 * log10(sum(clean.^2) / sum((clean - noise).^2));
fprintf('SNR after Bandpass Filtering (Fix #1): %.2f dB\n', SNR_after_bandpass);

%% LMS Filter
M = 12; 
mu_LMS = 0.1; 
w_LMS = zeros(M, 1); 
padded_signal = [zeros(M-1, 1); primary]; 
output_LMS = zeros(minLength, 1);

for n = 1:minLength
    u_vect = padded_signal(n:n+M-1); 
    e = noise(n) - w_LMS' * u_vect; 
    w_LMS = w_LMS + mu_LMS * e * u_vect; 
    output_LMS(n) = w_LMS' * u_vect;
end

filtered_signal_LMS = primary - output_LMS; 

%% NLMS Filter
mu_NLMS = 1; 
w_NLMS = zeros(M, 1); 
output_NLMS = zeros(minLength, 1);
Eps = 0.0001; 

for n = 1:minLength
    u_vect = padded_signal(n:n+M-1);
    mu_adapt = mu_NLMS / (Eps + norm(u_vect)^2);
    e = noise(n) - w_NLMS' * u_vect;
    w_NLMS = w_NLMS + mu_adapt * e * u_vect;
    output_NLMS(n) = w_NLMS' * u_vect;
end

filtered_signal_NLMS = primary - output_NLMS; 

%% RLS Filter
lambda = 1; % Forgetting factor  1 - 1 / (0.1 * M) =  0.1666 ; 
delta = 0.01; 
P = 1 / delta * eye(M); 
w_RLS = zeros(M, 1);
output_RLS = zeros(minLength, 1);

for n = 1:minLength
    u_vect = padded_signal(n:n+M-1); 
    PI = P * u_vect; 
    gain_k = PI / (lambda + u_vect' * PI); 
    e = noise(n) - w_RLS' * u_vect; 
    w_RLS = w_RLS + gain_k * e; 
    P = P / lambda - gain_k * (u_vect' * P) / lambda; 
    output_RLS(n) = w_RLS' * u_vect;
end

filtered_signal_RLS = primary - output_RLS; 

%% Calculate SNR Improvement
filtered_SNR_LMS = 10 * log10(sum(clean.^2) / sum((clean - filtered_signal_LMS).^2));
filtered_SNR_NLMS = 10 * log10(sum(clean.^2) / sum((clean - filtered_signal_NLMS).^2));
filtered_SNR_RLS = 10 * log10(sum(clean.^2) / sum((clean - filtered_signal_RLS).^2));

fprintf('SNR after LMS Filter: %.2f dB\n', mean(filtered_SNR_LMS));
fprintf('SNR after NLMS Filter: %.2f dB\n', mean(filtered_SNR_NLMS));
fprintf('SNR after RLS Filter: %.2f dB\n', mean(filtered_SNR_RLS));

%% Visualization of Results
figure;
subplot(5, 1, 1);
plot(clean);
title('Clean Signal');
xlim([0 minLength]);

subplot(5, 1, 2);
plot(primary);
title('Noisy Signal');
xlim([0 minLength]);

subplot(5, 1, 3);
plot(clean-filtered_signal_LMS);
title('Error signal (LMS)');
ylim([-1 1]);
xlim([0 minLength]);

subplot(5, 1, 4);
plot(clean-filtered_signal_NLMS);
title('Error signal (NLMS)');
ylim([-1 1]);
xlim([0 minLength]);

subplot(5, 1, 5);
plot(clean-filtered_signal_RLS);
title('Error signal (RLS)');
ylim([-1 1]);
xlim([0 minLength]);

tightfig();
saveas(gcf, 'Hospital_Ambient_Noises_NHS_1.pdf');

%% Spectrogram
windowLength = round(0.14*fs);  % 140 ms window length
hop_size = round(0.02*fs);      % 20 ms hop size
overlap = windowLength - hop_size; 
numBands = 128;                 % Number of Mel bands
hannWin = hann(windowLength, 'periodic'); 

% Compute Mel spectrogram for primary signal
[s_primary, f, t] = melSpectrogram(primary, fs, ...
                            "Window", kaiser(windowLength,18), ...
                            'OverlapLength', overlap, ...
                            'NumBands', numBands);

s_primary_db = 10 * log10(s_primary + eps);

% Compute Mel spectrogram for noise signal
[s_noise, ~, ~] = melSpectrogram(noise, fs, ...
                            "Window", kaiser(windowLength,18), ...
                            'OverlapLength', overlap, ...
                            'NumBands', numBands);

t = t - t(1); % Force time axis to start from 0

s_noise_db = 10 * log10(s_noise + eps);

% Plot Mel spectrograms
figure;

% Primary signal
subplot(2,1,1);
imagesc(t, f, s_primary_db);
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Primary Signal');
colorbar;
colormap jet;
set(gcf, 'Units', 'inches', 'Position', [1, 1, numel(primary)/fs, 1]);
ylim([0 4000]);
xlim([0 5]);

% Noise signal
subplot(2,1,2);
imagesc(t, f, s_noise_db);
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Noise Signal');
colorbar;
colormap jet;
set(gcf, 'Units', 'inches', 'Position', [1, 1, numel(noise)/fs, 1]);
ylim([0 4000]);
xlim([0 5]);

% Save figure
%tightfig();
%saveas(gcf, 'Mel_Spectrogram_real_heartbeatpdf');

%% Compare simulation synthethic heartbeat spectrogram with real data spectrogram

f1 = openfig('synthethic_heartbeat_melSpectrogram.fig', 'reuse'); % Load figure

% Get the axes from the loaded figure (this will give you only the axes, ignoring colorbars)
ax1 = findall(f1, 'Type', 'axes');  

length_of_spec_view = 10;

figure;
% --- Subplot 1: Synthetic Heartbeat Spectrogram ---
subplot(2,2,1);  % First subplot in a 2x2 grid
copyobj(allchild(ax1(2)), gca); 
title('Primary Signal (Synthetic Heartbeat)');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
xlim([0 length_of_spec_view]);
ylim([0 4000]);

% --- Subplot 2: Real Primary Signal Spectrogram ---
subplot(2,2,2);  % Second subplot in a 2x2 grid
imagesc(t, f, s_primary_db);   % Plot primary signal spectrogram
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Primary Signal (Real Heartbeat)');
xlim([0 length_of_spec_view]);
ylim([0 4000]);


% --- Subplot 3: Real Noise Signal Spectrogram (Copy) ---
subplot(2,2,3);  % Fourth subplot in a 2x2 grid
copyobj(allchild(ax1(1)), gca);  % Copy content from the second subplot
title('Noise Signal (Synthetic Heartbeat)');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
xlim([0 length_of_spec_view]);
ylim([0 4000]);

% --- Subplot 4: Real Noise Signal Spectrogram ---
subplot(2,2,4);  % Third subplot in a 2x2 grid
imagesc(t, f, s_noise_db);     % Plot noise signal spectrogram
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Noise Signal (Real Heartbeat)');
colormap jet;
xlim([0 length_of_spec_view]);
ylim([0 4000]);


% Set figure size
set(gcf, 'Units', 'inches', 'Position', [3.416666666666667,4.333333333333333,9.683333333333334,2.45]);


% Save the figure as PDF
tightfig();
saveas(gcf, 'Comparison_Mel_Spectrogram.pdf');


%% Cross-Correlation Before and After Noise Alignment
% figure;
% subplot(2,1,1);
% [xcorr_vals_before, lags_before] = xcorr(primary, noise);
% plot(lags_before/fs, xcorr_vals_before);
% title('Cross-Correlation Before Noise Alignment');
% xlabel('Time Lag (s)');
% ylabel('Cross-Correlation');
% 
% subplot(2,1,2);
% [xcorr_vals_after, lags_after] = xcorr(primary, noise);
% plot(lags_after/fs, xcorr_vals_after);
% title('Cross-Correlation After Noise Alignment');
% xlabel('Time Lag (s)');
% ylabel('Cross-Correlation');
% 
% saveas(gcf, 'Cross_Correlation_Alignment.pdf');

% %% Phase Difference Before and After Fix #3
% figure;
% subplot(2,1,1);
% plot(angle(analytic_primary ./ analytic_noise));
% title('Phase Difference Before Fix #3');
% xlabel('Sample Index');
% ylabel('Phase (radians)');
% 
% subplot(2,1,2);
% plot(phase_diff); % Now it shows the phase correction per sample
% title('Phase Difference After Fix #3');
% xlabel('Sample Index');
% ylabel('Phase (radians)');
% 
% saveas(gcf, 'Phase_Correction.pdf');

%% Export filtered signals
filtered_signal_LMS = filtered_signal_LMS / max(abs(filtered_signal_LMS));
filtered_signal_NLMS = filtered_signal_NLMS / max(abs(filtered_signal_NLMS));
filtered_signal_RLS = filtered_signal_RLS / max(abs(filtered_signal_RLS));

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

%% Save Figure trimmed
function tightfig()
    set(gcf, 'Units', 'Inches');
    pos = get(gcf, 'Position');
    set(gcf, 'PaperUnits', 'Inches');
    set(gcf, 'PaperSize', [pos(3) pos(4)]);
    set(gcf, 'PaperPosition', [0 0 pos(3) pos(4)]);
end
