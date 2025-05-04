clc;
clear;
close all;

%% Simulate Signal and Noise
[x, fs] = audioread('clean_synthetic_heartbeat.wav'); % Input signal
t = 0:1/fs:(length(x)-1)/fs; % Time vector
d = 0.5 * randn(size(t))'; % Simulated Gaussian noise
u = x + d; % Combine signal and noise

% Name of audio type to have in figures and folder name
suffix = 'Synthetic Heartbeat';

%% Define start time for trimming (1.5 seconds)
start_time = 5;  % in seconds
start_sample = round(start_time * fs);  % Convert time to sample index

%% Trim the signals to start from 1.5 seconds
u = u(start_sample:end);
d = d(start_sample:end);
x = x(start_sample:end);

% Find the minimum length
minLength = min([length(u), length(d), length(x)]);

% Calculate initial SNR
initial_SNR = 10 * log10(sum(x.^2) / sum((x - d).^2));  % The same as 10 * log10(sum(x.^2) / sum((x - d).^2))

initial_MSE = mean((x - d).^2);

%% Parameters
mu_values_LMS = [0.0001 0.001 0.002 0.005 0.0075 0.01 0.015 0.02 0.025 0.03 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
mu_values_NLMS = [0.0001 0.001 0.002 0.005 0.0075 0.01 0.015 0.02 0.025 0.03 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
lambda_values = [0.92 0.95 0.97 0.98 0.985 0.99 0.995 0.998 0.9985 0.999 0.9992 0.9995 0.9997 0.9999];
filter_length = [1 2 3 4 5 6 8 10 12 16 24 32 40 60 80 100]; 

totalIterations = length(filter_length) * (length(mu_values_LMS) + length(mu_values_NLMS)) + ...
                  length(filter_length) * length(lambda_values);

results(totalIterations) = struct('type', '', 'M', [], 'param', [], 'snr', -Inf);

parfor_progress(totalIterations); % Progress bar (if installed)

%% Main Loop
parfor idx = 1:totalIterations
    local_result = struct('type', '', 'M', [], 'param', [], 'snr', -Inf);

    totalLMS = length(filter_length) * length(mu_values_LMS);
    totalNLMS = length(filter_length) * length(mu_values_NLMS);
    totalMu = totalLMS + totalNLMS;

    if idx <= totalMu
        if idx <= totalLMS
            % LMS
            M_idx = ceil(idx / length(mu_values_LMS));
            mu_idx = mod(idx - 1, length(mu_values_LMS)) + 1;
            M = filter_length(M_idx);
            mu = mu_values_LMS(mu_idx);
            [y, e] = lms_filter(u, d, M, mu);
            filterName = 'LMS';
        else
            % NLMS
            adjusted_idx = idx - totalLMS;
            M_idx = ceil(adjusted_idx / length(mu_values_NLMS));
            mu_idx = mod(adjusted_idx - 1, length(mu_values_NLMS)) + 1;
            M = filter_length(M_idx);
            mu = mu_values_NLMS(mu_idx);
            [y, e] = nlms_filter(u, d, M, mu);
            filterName = 'NLMS';
        end

        if all(isfinite(y))
            snr_val = 10 * log10(sum(x.^2) / sum((x - e).^2));  % e = u - y
            local_result = struct('type', filterName, 'M', M, 'param', mu, 'snr', snr_val);
        end
    else
        % RLS
        idx_RLS = idx - totalMu;
        M_idx = ceil(idx_RLS / length(lambda_values));
        lambda_idx = mod(idx_RLS - 1, length(lambda_values)) + 1;

        M = filter_length(M_idx);
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
numOrders = length(filter_length);

bestSNR_LMS = -Inf(1, numOrders);
bestSNR_NLMS = -Inf(1, numOrders);
bestSNR_RLS = -Inf(1, numOrders);

optMu_LMS = zeros(1, numOrders);
optMu_NLMS = zeros(1, numOrders);
optLambda_RLS = zeros(1, numOrders);

for k = 1:totalIterations
    res = results(k);
    % 🛠 Map filter order (e.g., 4, 8...) to its index in filter_length array
    [isFound, M_idx] = ismember(res.M, filter_length);
    if  ~isFound
        warning("Filter order %d not found in filter_length!", res.M);
        continue;
    end

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
bestParamsLMS = zeros(length(filter_length), 3); % [Filter Length, Best SNR (LMS), Optimal mu (LMS)]
bestParamsNLMS = zeros(length(filter_length), 3); % [Filter Length, Best SNR (NLMS), Optimal mu (NLMS)]
bestParamsRLS = zeros(length(filter_length), 3); % [Filter Length, Best SNR (RLS), Optimal lambda (RLS)]

% Populate bestParamsLMS matrix
for M = 1:length(filter_length)
    bestParamsLMS(M, :) = [filter_length(M), bestSNR_LMS(M), optMu_LMS(M)];
end

% Populate bestParamsNLMS matrix
for M = 1:length(filter_length)
    bestParamsNLMS(M, :) = [filter_length(M), bestSNR_NLMS(M), optMu_NLMS(M)];
end

% Populate bestParamsRLS matrix
for M = 1:length(filter_length)
    bestParamsRLS(M, :) = [filter_length(M), bestSNR_RLS(M), optLambda_RLS(M)];
end

%% Display Best Parameters with Pretty Formatting

% Define header
header = 'Filter Length | Best SNR (LMS) | Optimal mu (LMS) | Optimal lambda (RLS)';

% Print the LMS parameters
fprintf('%-15s%-25s%-25s\n', 'Filter Length', 'Best SNR (LMS)', 'Optimal mu (LMS)');
for i = 1:length(filter_length)
    fprintf('%-15.4f%-25.4f%-25.4f\n', bestParamsLMS(i,1), bestParamsLMS(i,2), bestParamsLMS(i,3));
end

fprintf('\n'); % Newline for separation

% Print the NLMS parameters
fprintf('%-15s%-25s%-25s\n', 'Filter Length', 'Best SNR (NLMS)', 'Optimal mu (NLMS)');
for i = 1:length(filter_length)
    fprintf('%-15.4f%-25.4f%-25.4f\n', bestParamsNLMS(i,1), bestParamsNLMS(i,2), bestParamsNLMS(i,3));
end

fprintf('\n'); % Newline for separation

% Print the RLS parameters
fprintf('%-15s%-25s%-25s\n', 'Filter Length', 'Best SNR (RLS)', 'Optimal lambda (RLS)');
for i = 1:length(filter_length)
    fprintf('%-15.4f%-25.4f%-25.4f\n', bestParamsRLS(i,1), bestParamsRLS(i,2), bestParamsRLS(i,3));
end

%% Find Absolute Best Parameters for Each Algorithm

% LMS
[maxSNR_LMS, idx_LMS] = max(bestSNR_LMS);
bestFilterOrder_LMS = filter_length(idx_LMS);
bestMu_LMS = optMu_LMS(idx_LMS);

% NLMS
[maxSNR_NLMS, idx_NLMS] = max(bestSNR_NLMS);
bestFilterOrder_NLMS = filter_length(idx_NLMS);
bestMu_NLMS = optMu_NLMS(idx_NLMS);

% RLS
[maxSNR_RLS, idx_RLS] = max(bestSNR_RLS);
bestFilterOrder_RLS = filter_length(idx_RLS);
bestLambda_RLS = optLambda_RLS(idx_RLS);

%% Run the filters with the best parameter
[y_LMS, e_LMS, w_hist_LMS] = lms_filter(u, d, bestFilterOrder_LMS, bestMu_LMS);
[y_NLMS, e_NLMS, w_hist_NLMS] = nlms_filter(u, d, bestFilterOrder_NLMS, bestMu_NLMS);
[y_RLS, e_RLS, w_hist_RLS] = rls_filter(u, d, bestFilterOrder_RLS, bestLambda_RLS);

filtered_LMS = u - y_LMS;
filtered_NLMS = u - y_NLMS;
filtered_RLS = u - y_RLS;

% Calculate MSE value (a scalar)
MSE_LMS = mean((e_LMS).^2);
MSE_NLMS = mean((e_NLMS).^2);
MSE_RLS = mean((e_RLS).^2);

%% Display Results
fprintf('\n--- Absolute Best Parameters ---\n');

fprintf('Initial SNR: %.4f dB\n', initial_SNR);  % Print initial SNR
fprintf('Initial MSE: %.4f \n', initial_MSE);  % Print initial MSE

fprintf('\nLMS:\n');
fprintf('Best SNR: %.4f dB\n', maxSNR_LMS);
fprintf('SNR Improvment: %.4f dB\n', maxSNR_LMS-initial_SNR);
fprintf('Best MSE: %.4f \n', MSE_LMS);
fprintf('MSE Improvment: %.4f \n', initial_MSE-MSE_LMS);
fprintf('Filter Length: %d\n', bestFilterOrder_LMS);
fprintf('Optimal mu: %.4f\n', bestMu_LMS);

fprintf('\nNLMS:\n');
fprintf('Best SNR: %.4f dB\n', maxSNR_NLMS);
fprintf('SNR Improvment: %.4f dB\n', maxSNR_NLMS-initial_SNR);
fprintf('Best MSE: %.4f \n', MSE_NLMS);
fprintf('MSE Improvment: %.4f \n', initial_MSE-MSE_NLMS);
fprintf('Filter Length: %d\n', bestFilterOrder_NLMS);
fprintf('Optimal mu: %.4f\n', bestMu_NLMS);


fprintf('\nRLS:\n');
fprintf('Best SNR: %.4f dB\n', maxSNR_RLS);
fprintf('SNR Improvment: %.4f dB\n', maxSNR_RLS-initial_SNR);
fprintf('Best MSE: %.4f \n', MSE_RLS);
fprintf('MSE Improvment: %.4f \n', initial_MSE-MSE_RLS);
fprintf('Filter Length: %d\n', bestFilterOrder_RLS);
fprintf('Optimal lambda: %.4f\n', bestLambda_RLS);

%% Save simulation workspace

% Use current date for folder and filename
date_str = datestr(now, 'yyyy-dd-mm');  % e.g., '2025-19-04'

% Construct folder and filename
foldername = [date_str '_' suffix];
filename = fullfile(foldername, [foldername '.mat']);

% Create folder if it doesn't exist
if ~exist(foldername, 'dir')
    mkdir(foldername);
end

% % Save workspace inside that folder
% save(filename,'-v7.3');

%%%%% Save audio %%%%%

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
fprintf(fid,'Initial MSE: %.4f dB\n', initial_MSE);  % Print initial MSE

fprintf(fid, '\nLMS:\n');
fprintf(fid, 'Best SNR: %.4f dB\n', maxSNR_LMS);
fprintf(fid, 'Best MSE: %.4f dB\n', MSE_LMS);
fprintf(fid,'SNR Improvment: %.4f dB\n', maxSNR_LMS-initial_SNR);
fprintf(fid,'MSE Improvment: %.4f \n', initial_MSE-MSE_LMS);
fprintf(fid, 'Filter Length: %d\n', bestFilterOrder_LMS);
fprintf(fid, 'Optimal mu: %.4f\n', bestMu_LMS);

fprintf(fid, '\nNLMS:\n');
fprintf(fid, 'Best SNR: %.4f dB\n', maxSNR_NLMS);
fprintf(fid, 'Best MSE: %.4f dB\n', MSE_NLMS);
fprintf(fid,'SNR Improvment: %.4f dB\n', maxSNR_NLMS-initial_SNR);
fprintf(fid,'MSE Improvment: %.4f \n', initial_MSE-MSE_NLMS);
fprintf(fid, 'Filter Length: %d\n', bestFilterOrder_NLMS);
fprintf(fid, 'Optimal mu: %.4f\n', bestMu_NLMS);

fprintf(fid, '\nRLS:\n');
fprintf(fid, 'Best SNR: %.4f dB\n', maxSNR_RLS);
fprintf(fid, 'Best MSE: %.4f dB\n', MSE_RLS);
fprintf(fid,'SNR Improvment: %.4f dB\n', maxSNR_RLS-initial_SNR);
fprintf(fid,'MSE Improvment: %.4f \n', initial_MSE-MSE_RLS);
fprintf(fid, 'Filter Length: %d\n', bestFilterOrder_RLS);
fprintf(fid, 'Optimal lambda: %.4f\n', bestLambda_RLS);

fprintf(fid, '\n');

fprintf(fid, '---------------------------------\n');

% LMS Table
fprintf(fid, '%-15s%-25s%-25s\n', 'Filter Length', 'Best SNR (LMS)', 'Optimal mu (LMS)');
for i = 1:length(filter_length)
    fprintf(fid, '%-15.0f%-25.4f%-25.4f\n', bestParamsLMS(i,1), bestParamsLMS(i,2), bestParamsLMS(i,3));
end
fprintf(fid, '\n');

% NLMS Table
fprintf(fid, '%-15s%-25s%-25s\n', 'Filter Length', 'Best SNR (NLMS)', 'Optimal mu (NLMS)');
for i = 1:length(filter_length)
    fprintf(fid, '%-15.0f%-25.4f%-25.4f\n', bestParamsNLMS(i,1), bestParamsNLMS(i,2), bestParamsNLMS(i,3));
end
fprintf(fid, '\n');

% RLS Table
fprintf(fid, '%-15s%-25s%-25s\n', 'Filter Length', 'Best SNR (RLS)', 'Optimal lambda (RLS)');
for i = 1:length(filter_length)
    fprintf(fid, '%-15.0f%-25.4f%-25.4f\n', bestParamsRLS(i,1), bestParamsRLS(i,2), bestParamsRLS(i,3));
end

% Close file
fclose(fid);

%% Plot SNR vs Filter Length
figure;
plot(filter_length, bestSNR_LMS, '-o', 'LineWidth', 1.5, 'DisplayName', 'LMS');
hold on;
plot(filter_length, bestSNR_NLMS, '-s', 'LineWidth', 1.5, 'DisplayName', 'NLMS');
plot(filter_length, bestSNR_RLS, '-^', 'LineWidth', 1.5, 'DisplayName', 'RLS');
grid on;
xlabel('Filter Length (M)');
ylabel('Best Output SNR [dB]');
%title('Output SNR as a function of Filter Length');
legend('Location','best');

% Save figure
tightfig();
saveas(gcf, fullfile(foldername, ['Output SNR vs Filter Length - ' suffix '.pdf']));

%% Plot Mel Spectrogram for LMS, NLMS, and RLS in one figure

windowLength = round(0.14*fs);  % 140 ms window length
hop_size = round(0.01*fs);      % 10 ms hop size
overlap = windowLength - hop_size; 
numBands = 256;                 % Number of Mel bands
hannWin = hann(windowLength, 'periodic'); 

% Compute Mel spectrogram for clean signal
[s_x, f, t] = melSpectrogram(x, fs, ...
                            "Window", kaiser(windowLength,18), ...
                            'OverlapLength', overlap, ...
                            'NumBands', numBands);

s_x_db = 10 * log10(s_x + eps);



% Compute Mel spectrogram for primary signal (heart+noise)
[s_u, ~, ~] = melSpectrogram(u, fs, ...
                            "Window", kaiser(windowLength,18), ...
                            'OverlapLength', overlap, ...
                            'NumBands', numBands);

s_u_db = 10 * log10(s_u + eps);

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

% Find global min and max for color axis
cmin = -50;
cmax = max([max(s_x_db(:)), max(s_u_db(:)), max(s_filtered_LMS_db(:)), max(s_filtered_NLMS_db(:)), max(s_filtered_RLS_db(:))]);

% Plot Mel spectrograms in a single figure (5 subplots)
figure;

% Clean signal
subplot(5,1,1); 
imagesc(t, f, s_x_db); 
axis xy; 
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Clean Signal');
colorbar;
colormap jet;
caxis([cmin cmax]); % Set color axis
ylim([0 1000]);
xlim([0 5]);

% Primary signal (heart+noise)
subplot(5,1,2);
imagesc(t, f, s_u_db); 
axis xy; 
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Primary Signal (Clean + Noise)');
colorbar;
colormap jet;
caxis([cmin cmax]); % Set color axis
ylim([0 1000]);
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
caxis([cmin cmax]); % Set color axis
ylim([0 1000]);
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
caxis([cmin cmax]); % Set color axis
ylim([0 1000]);
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
caxis([cmin cmax]); % Set color axis
ylim([0 1000]);
xlim([0 5]);

% Set figure size and position
set(gcf, 'Units', 'inches', 'Position', [0, 0, 8.266666666666666, 9.866666666666667]);


% Save the figure
tightfig();

% Construct filename
saveas(gcf, fullfile(foldername, ['Mel Spectrogram - ' suffix '.pdf']));

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
title('Primary Signal (Clean + Noise)');
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
saveas(gcf, fullfile(foldername, ['Error signals - ' suffix '.pdf']));

% %% Frequency Response of Final Weights
% n_fft = 257;  % Number of FFT points
% x_range = linspace(0, fs/2, n_fft);  % Frequency axis from 0 to fs/2
% 
% figure;
% 
% % LMS
% subplot(3, 1, 1);
% w_LMS_final = w_hist_LMS(:, end);
% fr_LMS = 20 * log10(abs(freqz(w_LMS_final, 1, n_fft)));
% plot(x_range, fr_LMS);
% xlabel('Frequency (Hz)');
% ylabel('Magnitude (dB)');
% title('Frequency Response of LMS Filter');
% grid on;
% 
% % NLMS
% subplot(3, 1, 2);
% w_NLMS_final = w_hist_NLMS(:, end);
% fr_NLMS = 20 * log10(abs(freqz(w_NLMS_final, 1, n_fft)));
% plot(x_range, fr_NLMS);
% xlabel('Frequency (Hz)');
% ylabel('Magnitude (dB)');
% title('Frequency Response of NLMS Filter');
% grid on;
% 
% % RLS
% subplot(3, 1, 3);
% w_RLS_final = w_hist_RLS(:, end);
% fr_RLS = 20 * log10(abs(freqz(w_RLS_final, 1, n_fft)));
% plot(x_range, fr_RLS);
% xlabel('Frequency (Hz)');
% ylabel('Magnitude (dB)');
% title('Frequency Response of RLS Filter');
% grid on;
% 
% set(gcf, 'Units', 'inches', 'Position', [9.433333333333334,2.383333333333333,5.033333333333333,8.358333333333333]);
% 
% % Save figure
% tightfig();
% saveas(gcf, fullfile(foldername, ['Frequency Response of Final Weights - ' suffix '.pdf']));

% %% Create Learning curve
% 
% R = 100; % Number of runs for ensemble average
% N_plot = 1000; % Number of samples for MSE plotting
% mse_LMS_avg = zeros(N_plot, 1);  % To store average MSE for LMS
% mse_NLMS_avg = zeros(N_plot, 1); % To store average MSE for NLMS
% mse_RLS_avg = zeros(N_plot, 1);  % To store average MSE for RLS
% 
% for r = 1:R
%     % Run the filters
%     [y_LMS] = lms_filter(u, d, bestFilterOrder_LMS, bestMu_LMS);
%     [y_NLMS] = nlms_filter(u, d, bestFilterOrder_NLMS, bestMu_NLMS);
%     [y_RLS] = rls_filter(u, d, bestFilterOrder_RLS, bestLambda_RLS);
% 
%     % Calculate MSE (Mean Squared Error) at each iteration for each filter
%     mse_LMS = (u(1:N_plot)-y_LMS(1:N_plot)).^2;  % Squared error for LMS
%     mse_NLMS = (u(1:N_plot)-y_NLMS(1:N_plot)).^2;  % Squared error for NLMS
%     mse_RLS = (u(1:N_plot)-y_RLS(1:N_plot)).^2;  % Squared error for RLS
% 
%     % Accumulate MSE values for averaging
%     mse_LMS_avg = mse_LMS_avg + mse_LMS;
%     mse_NLMS_avg = mse_NLMS_avg + mse_NLMS;
%     mse_RLS_avg = mse_RLS_avg + mse_RLS;
% end
% 
% % Average MSE over all runs
% mse_LMS_avg = mse_LMS_avg / R;
% mse_NLMS_avg = mse_NLMS_avg / R;
% mse_RLS_avg = mse_RLS_avg / R;
% 
% %% Plot Learning Curve in 3 Subplots
% figure; set(gcf, 'Color', 'w');
% 
% % LMS Learning Curve
% subplot(3, 1, 1);  % Create subplot 1
% plot(mse_LMS_avg, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.2);
% xlabel('Iterations');  % Label as iterations since we're tracking filter updates
% ylabel('Avg MSE');
% legend('LMS', 'Location', 'northeast');
% %axis([0 1000 0 0.5])
% grid on;
% 
% % NLMS Learning Curve
% subplot(3, 1, 2);  % Create subplot 2
% plot(mse_NLMS_avg, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.2);
% xlabel('Iterations');  % Label as iterations
% ylabel('Avg MSE');
% legend('NLMS', 'Location', 'northeast');
% %axis([0 1000 0 0.5])
% grid on;
% 
% % RLS Learning Curve
% subplot(3, 1, 3);  % Create subplot 3
% plot(mse_RLS_avg, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 1.2);
% xlabel('Iterations');  % Label as iterations
% ylabel('Avg MSE');
% legend('RLS', 'Location', 'northeast');
% %axis([0 1000 0 0.5])
% grid on;
% 
% % Adjust layout
% tightfig();  % Ensure tight layout
% 
% % Save the learning curve plot
% saveas(gcf, fullfile(foldername, ['Learning Curve - ' suffix '.pdf']));

%% Error convergence curve comparison with smoothing

[y_LMS, e_LMS, w_hist_LMS] = lms_filter(u, d, bestFilterOrder_LMS, bestMu_LMS);
[y_NLMS, e_NLMS, w_hist_NLMS] = nlms_filter(u, d, bestFilterOrder_NLMS, bestMu_NLMS);
[y_RLS, e_RLS, w_hist_RLS] = rls_filter(u, d, bestFilterOrder_RLS, bestLambda_RLS);

N_plot = 2000; % Plotting only a portion
smooth_window = 50; % Moving average window size

% Calculate instantaneous squared error 
mse_inst_LMS = (d(1:N_plot)-y_LMS(1:N_plot)).^2;
mse_inst_NLMS = (d(1:N_plot)-y_NLMS(1:N_plot)).^2;
mse_inst_RLS = (d(1:N_plot)-y_RLS(1:N_plot)).^2;

% Smooth MSE using moving average
mse_curve_LMS_smooth = movmean(mse_inst_LMS, smooth_window);
mse_curve_NLMS_smooth = movmean(mse_inst_NLMS, smooth_window);
mse_curve_RLS_smooth = movmean(mse_inst_RLS, smooth_window);

% Find the maximum MSE across all filters
max_mse = max([max(mse_curve_LMS_smooth), max(mse_curve_NLMS_smooth), max(mse_curve_RLS_smooth)]);

% Plot Learning Curve in 3 Subplots
figure; set(gcf, 'Color', 'w');

% LMS Learning Curve
subplot(3, 1, 1);  % Create subplot 1
plot(mse_curve_LMS_smooth, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.2);
xlabel('Samples');
ylabel('MSE');
legend('LMS', 'Location', 'northeast');
grid on;
ylim([0 max_mse]);  % Set y-axis to have the same max limit

% NLMS Learning Curve
subplot(3, 1, 2);  % Create subplot 2
plot(mse_curve_NLMS_smooth, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.2);
xlabel('Samples');
ylabel('MSE');
legend('NLMS', 'Location', 'northeast');
grid on;
ylim([0 max_mse]);  % Set y-axis to have the same max limit

% RLS Learning Curve
subplot(3, 1, 3);  % Create subplot 3
plot(mse_curve_RLS_smooth, 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 1.2);
xlabel('Samples');
ylabel('MSE');
legend('RLS', 'Location', 'northeast');
grid on;
ylim([0 max_mse]);  % Set y-axis to have the same max limit

% Adjust layout
tightfig();  % Ensure tight layout

% Save the learning curve plot
saveas(gcf, fullfile(foldername, ['Error convergence curve - ' suffix '.pdf']));

%% LMS Filter Function
function [y, e, w_hist] = lms_filter(primary, secondary, M, mu)
    N = length(primary);
    w = zeros(M, 1);
    y = zeros(N,1);
    e = zeros(N,1);
    u_padded = [zeros(M-1, 1); secondary];
    w_hist = zeros(M, N); % Store filter weights

    for n = 1:N
        u_vec = u_padded(n:n+M-1);  % Current input vector
        y(n) = w' * u_vec;  % Filtered output
        e(n) = primary(n) - y(n); % Error signal
        w = w + mu * e(n) * u_vec; % Update weights

        w_hist(:, n) = w;
    end
end

%% NLMS Filter Function
function [y, e, w_hist] = nlms_filter(primary, secondary, M, mu)
    N = length(secondary);
    w = zeros(M, 1);
    y = zeros(N,1);
    e = zeros(N,1);
    eps = 0.1; % Stability constant 
    u_padded = [zeros(M-1, 1); secondary];
    w_hist = zeros(M, N); % Store filter weights

    for n = 1:N
        u_vec = u_padded(n:n+M-1); % Current input vector
        mu1 = mu / (eps + norm(u_vec)^2); % Normalized step size
        y(n) = w' * u_vec; % Filtered output
        e(n) = primary(n) - y(n); % Error signal
        w = w + mu1 * e(n) * u_vec; % Update weights

        w_hist(:, n) = w;
    end
end

%% RLS Filter Function
function [y, e, w_hist] = rls_filter(primary, secondary, M, lambda)
 N = length(secondary);
    w = zeros(M, 1); % Initialize filter weights
    y = zeros(N,1);
    e = zeros(N,1);
    delta = 0.01; % Initialization constant
    P = (1 / delta) * eye(M); % Initialize inverse correlation matrix
    u_padded = [sqrt(delta) * randn(M-1, 1); secondary]; % Pad noisy signal
    w_hist = zeros(M, N); % Store filter weights

    for n = 1:N
        u_vec = u_padded(n:n+M-1); % Current input vector
        PI = P * u_vec; % Intermediate calculation
        k = PI / (lambda + u_vec' * PI);  % Gain
        y(n) = w' * u_vec; % Filtered output
        e(n) = primary(n) - y(n); % Error signal
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