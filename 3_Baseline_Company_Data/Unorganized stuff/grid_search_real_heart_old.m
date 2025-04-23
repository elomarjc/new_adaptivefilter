clc;
clear;
close all;


%% Define the path for the data
[u, fs] = audioread("C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Data from AI heathway\Data_ANC\Experiment_Data\Hospital Ambient Noises\NHS\1\primary.wav");   %noise + clean signal
[d, ~] = audioread("C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Data from AI heathway\Data_ANC\Experiment_Data\Hospital Ambient Noises\NHS\1\secondary.wav");
[x, ~] = audioread("C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Data from AI heathway\Data_ANC\Experiment_Data\Hospital Ambient Noises\NHS\1\ZCH0019.wav");

%% Define start time for trimming (1.5 seconds)
start_time = 1.7;  % in seconds
start_sample = round(start_time * fs);  % Convert time to sample index

%% Trim the signals to start from 1.5 seconds
u = u(start_sample:end);
d = d(start_sample:end);
x = x(start_sample:end);

% Find the minimum length
minLength = min([length(u), length(d), length(x)]);

%% Bandpass Filtering to Remove Distortions
lowCutoff = 5; % Low cut frequency in Hz
highCutoff = 800; % High cut frequency in Hz
[b, a] = butter(4, [lowCutoff, highCutoff] / (fs/2), 'bandpass');
u = filtfilt(b, a, u);
d = filtfilt(b, a, d);

% Calculate initial SNR
initial_SNR = 10 * log10(sum(x.^2) / sum((x - d).^2));  % The same as 10 * log10(sum(x.^2) / sum((x - d).^2))

%% Parameters

% Rough range - tid 10.52
mu_values_LMS  = 0:0.05:1;
mu_values_NLMS = 0:0.05:1;
filterOrders = 1:20;
lambda_values  = 0:0.05:1;
enabled = 0;

% Fine range - tid 10.52
% mu_values_LMS  = 0.0005:0.001:0.03;
% mu_values_NLMS = 0.96:0.005:1;
% filterOrders = 1:25;
% lambda_values  = 0.0005:0.01:1;   %0.8:0.01:1;
%enabled = 1;

totalIterations = length(filterOrders) * (length(mu_values_LMS) + length(mu_values_NLMS)) + ...
                  length(filterOrders) * length(lambda_values);

results(totalIterations) = struct('type', '', 'M', [], 'param', [], 'snr', -Inf);

parfor_progress(totalIterations); % Progress bar (if installed)

%% Main Loop
parfor idx = 1:totalIterations
    local_result = struct('type', '', 'M', [], 'param', [], 'snr', -Inf);

    totalLMS = length(filterOrders) * length(mu_values_LMS);
    totalNLMS = length(filterOrders) * length(mu_values_NLMS);
    totalMu = totalLMS + totalNLMS;

    if idx <= totalMu
        if idx <= totalLMS
            % LMS
            M_idx = ceil(idx / length(mu_values_LMS));
            mu_idx = mod(idx - 1, length(mu_values_LMS)) + 1;
            M = filterOrders(M_idx);
            mu = mu_values_LMS(mu_idx);
            [y, e] = lms_filter(u, d, M, mu);
            filterName = 'LMS';
        else
            % NLMS
            adjusted_idx = idx - totalLMS;
            M_idx = ceil(adjusted_idx / length(mu_values_NLMS));
            mu_idx = mod(adjusted_idx - 1, length(mu_values_NLMS)) + 1;
            M = filterOrders(M_idx);
            mu = mu_values_NLMS(mu_idx);
            [y, e] = nlms_filter(u, d, M, mu);
            filterName = 'NLMS';
        end
        if all(isfinite(y))
            snr_val = 10 * log10(sum(x.^2) / sum((x - e).^2));
            local_result = struct('type', filterName, 'M', M, 'param', mu, 'snr', snr_val);
        end
    else
        % RLS
        idx_RLS = idx - totalMu;
        M_idx = ceil(idx_RLS / length(lambda_values));
        lambda_idx = mod(idx_RLS - 1, length(lambda_values)) + 1;

        M = filterOrders(M_idx);
        lambda_val = lambda_values(lambda_idx);

        [y, e] = rls_filter(u, d, M, lambda_val);

        if all(isfinite(y))
            snr_val = 10 * log10(sum(x.^2) / sum((x - e).^2));
            local_result = struct('type', 'RLS', 'M', M, 'param', lambda_val, 'snr', snr_val);
        end
    end

    results(idx) = local_result;
    parfor_progress;
end

parfor_progress(0);

%% Find Best Parameters
bestSNR_LMS = -Inf(1, length(filterOrders));
bestSNR_NLMS = -Inf(1, length(filterOrders));
bestSNR_RLS = -Inf(1, length(filterOrders));

optMu_LMS = zeros(1, length(filterOrders));
optMu_NLMS = zeros(1, length(filterOrders));
optLambda_RLS = zeros(1, length(filterOrders));

for k = 1:totalIterations
    res = results(k);
    M_idx = res.M;

    if strcmp(res.type, 'LMS')
        if res.snr > bestSNR_LMS(M_idx)
            bestSNR_LMS(M_idx) = res.snr;
            optMu_LMS(M_idx) = res.param;
        end
    elseif strcmp(res.type, 'NLMS')
        if res.snr > bestSNR_NLMS(M_idx)
            bestSNR_NLMS(M_idx) = res.snr;
            optMu_NLMS(M_idx) = res.param;
        end
    elseif strcmp(res.type, 'RLS')
        if res.snr > bestSNR_RLS(M_idx)
            bestSNR_RLS(M_idx) = res.snr;
            optLambda_RLS(M_idx) = res.param;
        end
    end
end

%% Create Matrix for Best Parameters
% Pre-allocate matrices for best parameters for each algorithm
bestParamsLMS = zeros(length(filterOrders), 3); % [Filter Order, Best SNR (LMS), Optimal mu (LMS)]
bestParamsNLMS = zeros(length(filterOrders), 3); % [Filter Order, Best SNR (NLMS), Optimal mu (NLMS)]
bestParamsRLS = zeros(length(filterOrders), 3); % [Filter Order, Best SNR (RLS), Optimal lambda (RLS)]

% Populate bestParamsLMS matrix
for M = 1:length(filterOrders)
    bestParamsLMS(M, :) = [filterOrders(M), bestSNR_LMS(M), optMu_LMS(M)];
end

% Populate bestParamsNLMS matrix
for M = 1:length(filterOrders)
    bestParamsNLMS(M, :) = [filterOrders(M), bestSNR_NLMS(M), optMu_NLMS(M)];
end

% Populate bestParamsRLS matrix
for M = 1:length(filterOrders)
    bestParamsRLS(M, :) = [filterOrders(M), bestSNR_RLS(M), optLambda_RLS(M)];
end

%% Display Best Parameters with Pretty Formatting

% Define header
header = 'Filter Order | Best SNR (LMS) | Optimal mu (LMS) | Optimal lambda (RLS)';

% Print the LMS parameters
fprintf('%-15s%-25s%-25s\n', 'Filter Order', 'Best SNR (LMS)', 'Optimal mu (LMS)');
for i = 1:length(filterOrders)
    fprintf('%-15.4f%-25.4f%-25.4f\n', bestParamsLMS(i,1), bestParamsLMS(i,2), bestParamsLMS(i,3));
end

fprintf('\n'); % Newline for separation

% Print the NLMS parameters
fprintf('%-15s%-25s%-25s\n', 'Filter Order', 'Best SNR (NLMS)', 'Optimal mu (NLMS)');
for i = 1:length(filterOrders)
    fprintf('%-15.4f%-25.4f%-25.4f\n', bestParamsNLMS(i,1), bestParamsNLMS(i,2), bestParamsNLMS(i,3));
end

fprintf('\n'); % Newline for separation

% Print the RLS parameters
fprintf('%-15s%-25s%-25s\n', 'Filter Order', 'Best SNR (RLS)', 'Optimal lambda (RLS)');
for i = 1:length(filterOrders)
    fprintf('%-15.4f%-25.4f%-25.4f\n', bestParamsRLS(i,1), bestParamsRLS(i,2), bestParamsRLS(i,3));
end

%% Find Absolute Best Parameters for Each Algorithm

% LMS
[maxSNR_LMS, idx_LMS] = max(bestSNR_LMS);
bestFilterOrder_LMS = filterOrders(idx_LMS);
bestMu_LMS = optMu_LMS(idx_LMS);

% NLMS
[maxSNR_NLMS, idx_NLMS] = max(bestSNR_NLMS);
bestFilterOrder_NLMS = filterOrders(idx_NLMS);
bestMu_NLMS = optMu_NLMS(idx_NLMS);

% RLS
[maxSNR_RLS, idx_RLS] = max(bestSNR_RLS);
bestFilterOrder_RLS = filterOrders(idx_RLS);
bestLambda_RLS = optLambda_RLS(idx_RLS);

%% Display Results
fprintf('\n--- Absolute Best Parameters ---\n');

fprintf('Initial SNR: %.4f dB\n', initial_SNR);  % Print initial SNR

fprintf('\nLMS:\n');
fprintf('Best SNR: %.4f dB\n', maxSNR_LMS);
fprintf('SNR Improvment: %.4f dB\n', maxSNR_LMS-initial_SNR);
fprintf('Filter Order: %d\n', bestFilterOrder_LMS);
fprintf('Optimal mu: %.4f\n', bestMu_LMS);

fprintf('\nNLMS:\n');
fprintf('Best SNR: %.4f dB\n', maxSNR_NLMS);
fprintf('SNR Improvment: %.4f dB\n', maxSNR_NLMS-initial_SNR);
fprintf('Filter Order: %d\n', bestFilterOrder_NLMS);
fprintf('Optimal mu: %.4f\n', bestMu_NLMS);

fprintf('\nRLS:\n');
fprintf('Best SNR: %.4f dB\n', maxSNR_RLS);
fprintf('SNR Improvment: %.4f dB\n', maxSNR_RLS-initial_SNR);
fprintf('Filter Order: %d\n', bestFilterOrder_RLS);
fprintf('Optimal lambda: %.4f\n', bestLambda_RLS);

%% Save simulation workspace

% Use current date for folder and filename
date_str = datestr(now, 'yyyy-dd-mm');  % e.g., '2025-19-04'

% Determine suffix based on enabled value
if enabled == 1
    suffix = 'fine';
else
    suffix = 'rough';
end

% Construct folder and filename
foldername = ['NoMSE_' date_str '_' suffix];
filename = fullfile(foldername, [foldername '.mat']);

% Create folder if it doesn't exist
if ~exist(foldername, 'dir')
    mkdir(foldername);
end

% % Save workspace inside that folder
% save(filename,'-v7.3');

%%%%% Save audio %%%%%
[y_LMS, ~, ~] = lms_filter(u, d, bestFilterOrder_LMS, bestMu_LMS);
[y_NLMS, ~, ~] = nlms_filter(u, d, bestFilterOrder_NLMS, bestMu_NLMS);
[y_RLS, ~, ~] = rls_filter(u, d, bestFilterOrder_RLS, bestLambda_RLS);

filtered_LMS = u - y_LMS;
filtered_NLMS = u - y_NLMS;
filtered_RLS = u - y_RLS;

% Normalize filtered signals (To prevent "Warning: Data clipped when writing file")
output_LMS_normalized = filtered_LMS / max(abs(filtered_LMS));
output_NLMS_normalized = filtered_NLMS / max(abs(filtered_NLMS));
output_RLS_normalized = filtered_RLS / max(abs(filtered_RLS));

% Save Filtered Signals in the newly created folder
audiowrite(fullfile(foldername, 'Filtered_LMS.wav'), output_LMS_normalized, fs);
audiowrite(fullfile(foldername, 'Filtered_NLMS.wav'), output_NLMS_normalized, fs);
audiowrite(fullfile(foldername, 'Filtered_RLS.wav'), output_RLS_normalized, fs);

fprintf('Filtered signals saved as audio files.\n\n');


%%%%%% Save parameters %%%%%

% Define result text file path (same folder as saved .mat)
results_txt_path = fullfile(foldername, 'best_parameters_summary.txt');

% Open file for writing (overwrite mode)
fid = fopen(results_txt_path, 'w');

% Check if file opened successfully
if fid == -1
    error('Failed to create result text file.');
end

% Summary Section
fprintf(fid, '--- Absolute Best Parameters ---\n');

fprintf(fid,'Initial SNR: %.4f dB\n', initial_SNR);  % Print initial SNR

fprintf(fid, '\nLMS:\n');
fprintf(fid, 'Best SNR: %.4f dB\n', maxSNR_LMS);
fprintf(fid,'SNR Improvment: %.4f dB\n', maxSNR_LMS-initial_SNR);
fprintf(fid, 'Filter Order: %d\n', bestFilterOrder_LMS);
fprintf(fid, 'Optimal mu: %.4f\n', bestMu_LMS);

fprintf(fid, '\nNLMS:\n');
fprintf(fid, 'Best SNR: %.4f dB\n', maxSNR_NLMS);
fprintf(fid,'SNR Improvment: %.4f dB\n', maxSNR_NLMS-initial_SNR);
fprintf(fid, 'Filter Order: %d\n', bestFilterOrder_NLMS);
fprintf(fid, 'Optimal mu: %.4f\n', bestMu_NLMS);

fprintf(fid, '\nRLS:\n');
fprintf(fid, 'Best SNR: %.4f dB\n', maxSNR_RLS);
fprintf(fid,'SNR Improvment: %.4f dB\n', maxSNR_RLS-initial_SNR);
fprintf(fid, 'Filter Order: %d\n', bestFilterOrder_RLS);
fprintf(fid, 'Optimal lambda: %.4f\n', bestLambda_RLS);

fprintf(fid, '\n');

fprintf(fid, '---------------------------------\n');

% LMS Table
fprintf(fid, '%-15s%-25s%-25s\n', 'Filter Order', 'Best SNR (LMS)', 'Optimal mu (LMS)');
for i = 1:length(filterOrders)
    fprintf(fid, '%-15.0f%-25.4f%-25.4f\n', bestParamsLMS(i,1), bestParamsLMS(i,2), bestParamsLMS(i,3));
end
fprintf(fid, '\n');

% NLMS Table
fprintf(fid, '%-15s%-25s%-25s\n', 'Filter Order', 'Best SNR (NLMS)', 'Optimal mu (NLMS)');
for i = 1:length(filterOrders)
    fprintf(fid, '%-15.0f%-25.4f%-25.4f\n', bestParamsNLMS(i,1), bestParamsNLMS(i,2), bestParamsNLMS(i,3));
end
fprintf(fid, '\n');

% RLS Table
fprintf(fid, '%-15s%-25s%-25s\n', 'Filter Order', 'Best SNR (RLS)', 'Optimal lambda (RLS)');
for i = 1:length(filterOrders)
    fprintf(fid, '%-15.0f%-25.4f%-25.4f\n', bestParamsRLS(i,1), bestParamsRLS(i,2), bestParamsRLS(i,3));
end

% Close file
fclose(fid);

%% Plot SNR vs Filter Order
figure;
plot(filterOrders, bestSNR_LMS, '-o', 'LineWidth', 1.5, 'DisplayName', 'LMS');
hold on;
plot(filterOrders, bestSNR_NLMS, '-s', 'LineWidth', 1.5, 'DisplayName', 'NLMS');
plot(filterOrders, bestSNR_RLS, '-^', 'LineWidth', 1.5, 'DisplayName', 'RLS');
grid on;
xlabel('Filter Order (M)');
ylabel('Best Output SNR [dB]');
%title('Output SNR as a function of Filter Length');
legend('Location','best');

% Save figure
tightfig();
saveas(gcf, fullfile(foldername, 'Output SNR vs Filter Length - Real Heart.pdf'));

%% Plot Mel Spectrogram for LMS, NLMS, and RLS in one figure

windowLength = round(0.14*fs);  % 140 ms window length
hop_size = round(0.02*fs);      % 20 ms hop size
overlap = windowLength - hop_size; 
numBands = 128;                 % Number of Mel bands
hannWin = hann(windowLength, 'periodic'); 

% Compute Mel spectrogram for primary signal
[s_u, f, t] = melSpectrogram(u, fs, ...
                            "Window", kaiser(windowLength,18), ...
                            'OverlapLength', overlap, ...
                            'NumBands', numBands);

s_u_db = 10 * log10(s_u + eps);

% Compute Mel spectrogram for noise signal
[s_d, ~, ~] = melSpectrogram(d, fs, ...
                            "Window", kaiser(windowLength,18), ...
                            'OverlapLength', overlap, ...
                            'NumBands', numBands);

s_d_db = 10 * log10(s_d + eps);

t = t - t(1); % Force time axis to start from 0

% Compute Mel spectrogram for filtered LMS signal
[s_filtered_LMS, ~, ~] = melSpectrogram(filtered_LMS, fs, ...
                            "Window", kaiser(windowLength,18), ...
                            'OverlapLength', overlap, ...
                            'NumBands', numBands);

s_filtered_LMS_db = 10 * log10(s_filtered_LMS + eps);

% Compute Mel spectrogram for filtered NLMS signal
[s_filtered_NLMS, ~, ~] = melSpectrogram(filtered_NLMS, fs, ...
                            "Window", kaiser(windowLength,18), ...
                            'OverlapLength', overlap, ...
                            'NumBands', numBands);

s_filtered_NLMS_db = 10 * log10(s_filtered_NLMS + eps);

% Compute Mel spectrogram for filtered RLS signal
[s_filtered_RLS, ~, ~] = melSpectrogram(filtered_RLS, fs, ...
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


% Save the figure
tightfig();
saveas(gcf, fullfile(foldername, 'Mel_Spectrogram - Real Heart.pdf'));

%% Plot Error signals

figure;
subplot(5, 1, 1);
plot(x);
title('Clean Signal');
xlabel('Sample Number');
ylabel('Amplitude');
ylim([-1.1 1.1]);
xlim([0 length(u)]);

subplot(5, 1, 2);
plot(u);
title('Noisy Signal');
xlabel('Sample Number');
ylabel('Amplitude');
ylim([-1.1 1.1]);
xlim([0 length(u)]);

subplot(5, 1, 3); 
plot(x-filtered_LMS);
title('Error signal (LMS)');
xlabel('Sample Number');
ylabel('Amplitude');
ylim([-1.1 1.1]);
xlim([0 length(u)]);

subplot(5, 1, 4);
plot(x-filtered_NLMS);
title('Error signal (NLMS)');
xlabel('Sample Number');
ylabel('Amplitude');
ylim([-1.1 1.1]);
xlim([0 length(u)]);

subplot(5, 1, 5);
plot(x-filtered_RLS);
title('Error signal (RLS)');
xlabel('Sample Number');
ylabel('Amplitude');
ylim([-1.1 1.1]);
xlim([0 length(u)]);

% Save figure
tightfig();
saveas(gcf, fullfile(foldername, 'Error signals - Real Heart.pdf'));

%% Frequency Response of Final Weights
n_fft = 257;  % Number of FFT points
x_range = linspace(0, fs/2, n_fft);  % Frequency axis from 0 to fs/2

[y_LMS, ~, w_hist_LMS] = lms_filter(u, d, bestFilterOrder_LMS, bestMu_LMS);
[y_NLMS, ~, w_hist_NLMS] = nlms_filter(u, d, bestFilterOrder_NLMS, bestMu_NLMS);
[y_RLS, ~, w_hist_RLS] = rls_filter(u, d, bestFilterOrder_RLS, bestLambda_RLS);

figure;

% LMS
subplot(3, 1, 1);
w_LMS_final = w_hist_LMS(:, end);
fr_LMS = 20 * log10(abs(freqz(w_LMS_final, 1, n_fft)));
plot(x_range, fr_LMS);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Frequency Response of LMS Filter');
grid on;

% NLMS
subplot(3, 1, 2);
w_NLMS_final = w_hist_NLMS(:, end);
fr_NLMS = 20 * log10(abs(freqz(w_NLMS_final, 1, n_fft)));
plot(x_range, fr_NLMS);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Frequency Response of NLMS Filter');
grid on;

% RLS
subplot(3, 1, 3);
w_RLS_final = w_hist_RLS(:, end);
fr_RLS = 20 * log10(abs(freqz(w_RLS_final, 1, n_fft)));
plot(x_range, fr_RLS);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Frequency Response of RLS Filter');
grid on;

set(gcf, 'Units', 'inches', 'Position', [9.433333333333334,2.383333333333333,5.033333333333333,8.358333333333333]);

% Save figure
tightfig();
saveas(gcf, fullfile(foldername, 'Frequency Response of Final Weights - Real Heart.pdf'));

%% LMS Filter Function
function [y, e, w_hist] = lms_filter(u, d, M, mu)
    N = length(u);
    w = zeros(M,1);
    y = zeros(N,1);
    e = zeros(N,1);
    u_padded = [zeros(M-1, 1); u];
    w_hist = zeros(M, N); % Store filter weights

    for n = 1:N
        u_vec = u_padded(n:n+M-1);  % Current input vector
        y(n) = w' * u_vec;  % Filtered output
        e(n) = d(n) - y(n); % Error signal
        w = w + mu * e(n) * u_vec; % Update weights
        
        w_hist(:, n) = w;
    end
end

%% NLMS Filter Function
function [y, e, w_hist] = nlms_filter(u, d, M, mu)
    N = length(u);
    w = zeros(M,1);
    y = zeros(N,1);
    e = zeros(N,1);
    eps = 0.0001; % Stability constant 
    u_padded = [zeros(M-1, 1); u];
    w_hist = zeros(M, N); % Store filter weights
    
    for n = 1:N
        u_vec = u_padded(n:n+M-1); % Current input vector
        y(n) = w' * u_vec; % Filtered output
        e(n) = d(n) - y(n); % Error signal
        mu1 = mu / (eps + norm(u_vec)^2); % Normalized step size
        w = w + mu1 * e(n) * u_vec; % Update weights
        
        w_hist(:, n) = w;
    end
end

%% RLS Filter Function
function [y, e, w_hist] = rls_filter(u, d, M, lambda)
    N = length(u);
    w = zeros(M,1); % Initialize filter weights
    y = zeros(N,1);
    e = zeros(N,1);
    delta = 0.01; % Initialization constant
    P = (1/delta) * eye(M); % Initialize inverse correlation matrix
    u_padded = [sqrt(delta) * randn(M-1, 1); u]; % Pad noisy signal
    w_hist = zeros(M, N); % Store filter weights
    
    for n = 1:N
        u_vec = u_padded(n:n+M-1); % Current input vector
        PI = P * u_vec; % Intermediate calculation
        k = PI / (lambda + u_vec' * PI);  % Gain
        y(n) = w' * u_vec; % Filtered output
        e(n) = d(n) - y(n); % Error signal
        w = w + k * e(n); % Update weights
        P = (P - k * u_vec' * P) / lambda; % Update inverse correlation matrix

        w_hist(:, n) = w;
    end
end

%% Save Figure trimmed
function tightfig()
    set(gcf, 'Units', 'Inches');
    pos = get(gcf, 'Position');
    set(gcf, 'PaperUnits', 'Inches');
    set(gcf, 'PaperSize', [pos(3) pos(4)]);
    set(gcf, 'PaperPosition', [0 0 pos(3) pos(4)]);
end