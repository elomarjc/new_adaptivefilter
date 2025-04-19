clc;
clear;
close all;

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

%% Parameters
mu_values_LMS = 0.001:0.1:1;   % Separate step size range for LMS
mu_values_NLMS = 0.001:0.1:1;    % Separate step size range for NLMS
filterOrders = 1:15;
lambda_values = 0.001:0.1:1;

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
            [y, ~] = lms_filter(noise, primary, M, mu);
            filterName = 'LMS';
        else
            % NLMS
            adjusted_idx = idx - totalLMS;
            M_idx = ceil(adjusted_idx / length(mu_values_NLMS));
            mu_idx = mod(adjusted_idx - 1, length(mu_values_NLMS)) + 1;
            M = filterOrders(M_idx);
            mu = mu_values_NLMS(mu_idx);
            [y, ~] = nlms_filter(noise, primary, M, mu);
            filterName = 'NLMS';
        end

        if all(isfinite(y))
            snr_val = snr(clean, clean - (primary - y));
            local_result = struct('type', filterName, 'M', M, 'param', mu, 'snr', snr_val);
        end
    else
        % RLS
        idx_RLS = idx - totalMu;
        M_idx = ceil(idx_RLS / length(lambda_values));
        lambda_idx = mod(idx_RLS - 1, length(lambda_values)) + 1;

        M = filterOrders(M_idx);
        lambda_val = lambda_values(lambda_idx);

        [y, ~] = rls_filter(noise, primary, M, lambda_val);

        if all(isfinite(y))
            snr_val = snr(clean, clean - (primary - y));
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

fprintf('\nLMS:\n');
fprintf('Best SNR: %.4f dB\n', maxSNR_LMS);
fprintf('Filter Order: %d\n', bestFilterOrder_LMS);
fprintf('Optimal mu: %.4f\n', bestMu_LMS);

fprintf('\nNLMS:\n');
fprintf('Best SNR: %.4f dB\n', maxSNR_NLMS);
fprintf('Filter Order: %d\n', bestFilterOrder_NLMS);
fprintf('Optimal mu: %.4f\n', bestMu_NLMS);

fprintf('\nRLS:\n');
fprintf('Best SNR: %.4f dB\n', maxSNR_RLS);
fprintf('Filter Order: %d\n', bestFilterOrder_RLS);
fprintf('Optimal lambda: %.4f\n', bestLambda_RLS);


%% Plot SNR vs Filter Order
figure;
plot(filterOrders, bestSNR_LMS, '-o', 'LineWidth', 1.5, 'DisplayName', 'LMS');
hold on;
plot(filterOrders, bestSNR_NLMS, '-s', 'LineWidth', 1.5, 'DisplayName', 'NLMS');
plot(filterOrders, bestSNR_RLS, '-^', 'LineWidth', 1.5, 'DisplayName', 'RLS');
grid on;
xlabel('Filter Order (M)');
ylabel('Best Output SNR [dB]');
title('Output SNR as a function of Filter Order');
legend('Location','northwest');

%% Plot Optimal mu / lambda vs Filter Order
figure;
plot(filterOrders, optMu_LMS, '-o', 'LineWidth', 1.5, 'DisplayName', 'Optimal \mu LMS');
hold on;
plot(filterOrders, optMu_NLMS, '-s', 'LineWidth', 1.5, 'DisplayName', 'Optimal \mu NLMS');
plot(filterOrders, optLambda_RLS, '-^', 'LineWidth', 1.5, 'DisplayName', 'Optimal \lambda RLS');
grid on;
xlabel('Filter Order (M)');
ylabel('Optimal Step Size / Forgetting Factor');
title('Optimal \mu / \lambda vs Filter Order');
legend('Location','best');

%% Plot Maximum SNR Achieved by Each Algorithm
figure;
maxSNRs = [max(bestSNR_LMS), max(bestSNR_NLMS), max(bestSNR_RLS)];
bar(maxSNRs);
set(gca, 'XTickLabel', {'LMS', 'NLMS', 'RLS'});
ylabel('Best Output SNR (dB)');
title('Maximum SNR Achieved by Each Algorithm');
grid on;

%% Plot Spectrogram
% windowLength = round(0.14*fs);  % 140 ms window length
% hop_size = round(0.02*fs);      % 20 ms hop size
% overlap = windowLength - hop_size; 
% numBands = 128;                 % Number of Mel bands
% hannWin = hann(windowLength, 'periodic'); 
% 
% % Compute Mel spectrogram for primary signal
% [s_primary, f, t] = melSpectrogram(primary, fs, ...
%                             "Window", kaiser(windowLength,18), ...
%                             'OverlapLength', overlap, ...
%                             'NumBands', numBands);
% 
% s_primary_db = 10 * log10(s_primary + eps);
% 
% % Compute Mel spectrogram for noise signal
% [s_noise, ~, ~] = melSpectrogram(noise, fs, ...
%                             "Window", kaiser(windowLength,18), ...
%                             'OverlapLength', overlap, ...
%                             'NumBands', numBands);
% 
% t = t - t(1); % Force time axis to start from 0
% 
% s_noise_db = 10 * log10(s_noise + eps);
% 
% % Plot Mel spectrograms
% figure;
% 
% % Primary signal
% subplot(2,1,1);
% imagesc(t, f, s_primary_db);
% axis xy;
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% title('Primary Signal');
% colorbar;
% colormap jet;
% set(gcf, 'Units', 'inches', 'Position', [3.416666666666667,4.333333333333333,9.683333333333334,2.4]);
% ylim([0 1000]);
% xlim([0 5]);
% 
% % Noise signal
% subplot(2,1,2);
% imagesc(t, f, s_noise_db);
% axis xy;
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% title('Noise Signal');
% colorbar;
% colormap jet;
% set(gcf, 'Units', 'inches', 'Position', [3.416666666666667,4.333333333333333,9.683333333333334,2.4]);
% ylim([0 1000]);
% xlim([0 5]);
% 
% % Save figure
% tightfig();
% saveas(gcf, 'Mel_Spectrogram_sinus_wave.pdf');

%%  Plot LMS Time evolution of the filter weights for the best SNR configuration

% Apply LMS filter with best parameters (we assume noise, primary, bestFilterOrder_LMS, bestMu_LMS are already defined)
[y_LMS, ~, w_hist_LMS] = lms_filter(noise, primary, bestFilterOrder_LMS, bestMu_LMS);

% Time steps to plot (choose a subset if needed to prevent clutter)
num_cycles = size(w_hist_LMS, 2); % Total number of cycles
subset_cycles = min(num_cycles, 500); % Display up to 500 cycles or adjust as needed

% Define a custom color map with more colors (we can use the 'hsv' colormap)
colors = hsv(bestFilterOrder_LMS); % 'hsv' generates a larger number of distinct colors

% Plot the evolution of filter weights for each tap
figure;
hold on;
for i = 1:bestFilterOrder_LMS
    plot(1:subset_cycles, w_hist_LMS(i, 1:subset_cycles), 'LineWidth', 1.5, 'Color', colors(i, :));
end
hold off;

grid on;
title('LMS Filter Weight Evolution (First 500 Cycles)');
xlabel('Cycle (n)');
ylabel('Weight Value');

% Generate the legend with distinct labels for each weight
legend_entries_lms = arrayfun(@(i) sprintf('w_{%d}(n)', i-1), 1:bestFilterOrder_LMS, 'UniformOutput', false);

% Adjust the font size and organize the legend in 3 columns
legend_handle__lms = legend(legend_entries_lms, 'Location', 'Best', 'FontSize', 8, 'NumColumns', 3);

%% Plot NLMS time evolution of the filter weights for the best SNR configuration
[y_NLMS, ~, w_hist_NLMS] = nlms_filter(noise, primary, bestFilterOrder_NLMS, bestMu_NLMS);

% Time steps to plot (choose a subset if needed to prevent clutter)
num_cycles = size(w_hist_NLMS, 2); % Total number of cycles
subset_cycles = min(num_cycles, 500); % Display up to 500 cycles or adjust as needed

% Define a custom color map with more colors (we can use the 'hsv' colormap)
colors = hsv(bestFilterOrder_NLMS); % 'hsv' generates a larger number of distinct colors

% Plot the evolution of filter weights for each tap
figure;
hold on;
for i = 1:bestFilterOrder_NLMS
    plot(1:subset_cycles, w_hist_NLMS(i, 1:subset_cycles), 'LineWidth', 1.5, 'Color', colors(i, :));
end
hold off;

grid on;
title('NLMS Filter Weight Evolution (First 500 Cycles)');
xlabel('Cycle (n)');
ylabel('Weight Value');

% Generate the legend with distinct labels for each weight
legend_entries_nlms = arrayfun(@(i) sprintf('w_{%d}(n)', i-1), 1:bestFilterOrder_NLMS, 'UniformOutput', false);

% Adjust the font size and organize the legend in 3 columns
legend_handle_nlms = legend(legend_entries_nlms, 'Location', 'Best', 'FontSize', 8, 'NumColumns', 3);

%% Plot RLS time evolution of the filter weights for the best SNR configuration
[y_RLS, ~, w_hist_RLS] = rls_filter(noise, primary, bestFilterOrder_RLS, bestLambda_RLS);

% Time steps to plot (choose a subset if needed to prevent clutter)
num_cycles = size(w_hist_RLS, 2); % Total number of cycles
subset_cycles = min(num_cycles, 500); % Display up to 500 cycles or adjust as needed

% Define a custom color map with more colors (we can use the 'hsv' colormap)
colors = hsv(bestFilterOrder_RLS); % 'hsv' generates a larger number of distinct colors

% Plot the evolution of filter weights for each tap
figure;
hold on;
for i = 1:bestFilterOrder_RLS
    plot(1:subset_cycles, w_hist_RLS(i, 1:subset_cycles), 'LineWidth', 1.5, 'Color', colors(i, :));
end
hold off;

grid on;
title('RLS Filter Weight Evolution (First 500 Cycles)');
xlabel('Cycle (n)');
ylabel('Weight Value');

% Generate the legend with distinct labels for each weight
legend_entries_rls = arrayfun(@(i) sprintf('w_{%d}(n)', i-1), 1:bestFilterOrder_RLS, 'UniformOutput', false);

% Adjust the font size and organize the legend in 3 columns
legend_handle_rls = legend(legend_entries_rls, 'Location', 'Best', 'FontSize', 8, 'NumColumns', 3);

%% Plot Error signals
figure;
subplot(5, 1, 1);
plot(clean);
title('Clean Signal');
xlabel('Sample Number');
ylabel('Amplitude');
ylim([-1.1 1.1]);

subplot(5, 1, 2);
plot(primary);
title('Noisy Signal');
xlabel('Sample Number');
ylabel('Amplitude');
ylim([-1.1 1.1]);

subplot(5, 1, 3); 
plot(clean-(primary-y_LMS));
title('Error signal (LMS)');
xlabel('Sample Number');
ylabel('Amplitude');
ylim([-1.1 1.1]);

subplot(5, 1, 4);
plot(clean-(primary-y_NLMS));
title('Error signal (NLMS)');
xlabel('Sample Number');
ylabel('Amplitude');
ylim([-1.1 1.1]);

subplot(5, 1, 5);
plot(clean-(primary-y_RLS));
title('Error signal (RLS)');
xlabel('Sample Number');
ylabel('Amplitude');
ylim([-1.1 1.1]);

% Save figure
tightfig();
saveas(gcf, 'sinus_wave_error_signals.pdf');

%% Frequency Response of Final Weights
n_fft = 257;  % Number of FFT points
x_range = linspace(0, fs/2, n_fft);  % Frequency axis from 0 to fs/2

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

%% LMS Filter Function
function [y, e, w_hist] = lms_filter(d, x, M, mu)
    N = length(x);
    w = zeros(M,1);
    y = zeros(N,1);
    e = zeros(N,1);
    x_padded = [zeros(M-1, 1); x];
    w_hist = zeros(M, N); % Store filter weights

    for n = 1:N
        x_vec = x_padded(n:n+M-1);  % Current input vector
        e = d(n) - w' * x_vec; % Error signal
        w = w + mu * e * x_vec; % Update weights
        y(n) = w' * x_vec;  % Filtered output
        
        w_hist(:, n) = w;
    end
end

%% NLMS Filter Function
function [y, e, w_hist] = nlms_filter(d, x, M, mu)
    N = length(x);
    w = zeros(M,1);
    y = zeros(N,1);
    e = zeros(N,1);
    eps = 0.0001; % Stability constant 
    x_padded = [zeros(M-1, 1); x];
    w_hist = zeros(M, N); % Store filter weights
    for n = 1:N
        x_vec = x_padded(n:n+M-1); % Current input vector
        mu1 = mu / (eps + norm(x_vec)^2); % Normalized step size
        e = d(n) - w' * x_vec; % Error signal
        w = w + mu1 * e * x_vec; % Update weights
        y(n) = w' * x_vec; % Filtered output
        
        w_hist(:, n) = w;
    end
end

%% RLS Filter Function
function [y, e, w_hist] = rls_filter(d, x, M, lambda)
    N = length(x);
    w = zeros(M,1); % Initialize filter weights
    y = zeros(N,1);
    e = zeros(N,1);
    delta = 0.01; % Initialization constant
    P = (1/delta) * eye(M); % Initialize inverse correlation matrix
    x_padded = [sqrt(delta) * randn(M-1, 1); x]; % Pad noisy signal
    w_hist = zeros(M, N); % Store filter weights
    for n = 1:N
        x_vec = x_padded(n:n+M-1); % Current input vector
        PI = P * x_vec; % Intermediate calculation
        k = PI / (lambda + x_vec' * PI);  % Gain
        e = d(n) - w' * x_vec; % Error signal
        w = w + e * k; % Update weights
        P = P / lambda - k * (x_vec' * P) / lambda; % Update inverse correlation matrix
        y(n) = w' * x_vec; % Filtered output
        
        w_hist(:, n) = w;
    end

end

%% Save simulation workspace

mu_values = mu_values_LMS;

% Extract step size info
mu_min  = min(mu_values);
mu_max  = max(mu_values);
mu_step = mu_values(2) - mu_values(1);

% Function to convert to string and remove trailing zeros
mu_to_str = @(x) regexprep(num2str(x, '%.10g'), '\.?0+$', '');

% Convert values to strings
mu_min_str  = strrep(mu_to_str(mu_min), '.', '');
mu_step_str = strrep(mu_to_str(mu_step), '.', '');
mu_max_str  = strrep(mu_to_str(mu_max), '.', '');

% Filter order bounds
order_min = min(filterOrders);
order_max = max(filterOrders);

% Final filename and folder name
filename_base = sprintf('FilterOrder_%d-%d_StepSizeRange_%s-%s-%s', ...
    order_min, order_max, mu_step_str, mu_min_str, mu_max_str);

foldername = filename_base;
filename = fullfile(foldername, [filename_base '.mat']);

% Create folder if it doesn't exist
if ~exist(foldername, 'dir')
    mkdir(foldername);
end

% Save workspace inside that folder
save(filename);

%% Save parameters

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

fprintf(fid, '\nLMS:\n');
fprintf(fid, 'Best SNR: %.4f dB\n', maxSNR_LMS);
fprintf(fid, 'Filter Order: %d\n', bestFilterOrder_LMS);
fprintf(fid, 'Optimal mu: %.4f\n', bestMu_LMS);

fprintf(fid, '\nNLMS:\n');
fprintf(fid, 'Best SNR: %.4f dB\n', maxSNR_NLMS);
fprintf(fid, 'Filter Order: %d\n', bestFilterOrder_NLMS);
fprintf(fid, 'Optimal mu: %.4f\n', bestMu_NLMS);

fprintf(fid, '\nRLS:\n');
fprintf(fid, 'Best SNR: %.4f dB\n', maxSNR_RLS);
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

%%
%%%%% Save audio %%%%%
[y_LMS, ~, ~] = lms_filter(noise, primary, bestFilterOrder_LMS, bestMu_LMS);
[y_NLMS, ~, ~] = nlms_filter(noise, primary, bestFilterOrder_NLMS, bestMu_NLMS);
[y_RLS, ~, ~] = rls_filter(noise, primary, bestFilterOrder_RLS, bestLambda_RLS);

filtered_LMS = (primary - y_LMS);
filtered_NLMS = (primary - y_NLMS);
filtered_RLS = (primary - y_RLS);

% Normalize filtered signals (To prevent "Warning: Data clipped when writing file")
output_LMS_normalized = (primary - y_LMS) / max(abs(filtered_LMS));
output_NLMS_normalized = (primary - y_NLMS) / max(abs(filtered_NLMS));
output_RLS_normalized = (primary - y_RLS) / max(abs(filtered_RLS));

% Save Filtered Signals in the newly created folder
audiowrite(fullfile(foldername, 'Filtered_LMS.wav'), output_LMS_normalized, fs);
audiowrite(fullfile(foldername, 'Filtered_NLMS.wav'), output_NLMS_normalized, fs);
audiowrite(fullfile(foldername, 'Filtered_RLS.wav'), output_RLS_normalized, fs);

fprintf('Filtered signals saved as audio files.\n\n');


%% Save Figure trimmed
function tightfig()
    set(gcf, 'Units', 'Inches');
    pos = get(gcf, 'Position');
    set(gcf, 'PaperUnits', 'Inches');
    set(gcf, 'PaperSize', [pos(3) pos(4)]);
    set(gcf, 'PaperPosition', [0 0 pos(3) pos(4)]);
end

%% LMS: Plot time evolution of the filter weights for the best SNR configuration
[y_LMS, ~, w_hist_LMS] = lms_filter(noise, primary, bestFilterOrder_LMS, bestMu_LMS);

% Plot weights evolution
figure;
plot(w_hist_LMS');
grid on;
title('LMS Weights Evolution');
xlabel('Cycle (n)');
ylabel('Magnitude');
legend(arrayfun(@(i) sprintf('w_%d(n)', i - 1), 1:bestFilterOrder_LMS, 'UniformOutput', false));

%% NLMS: Plot time evolution of the filter weights for the best SNR configuration
[y_NLMS, ~, w_hist_NLMS] = nlms_filter(noise, primary, bestFilterOrder_NLMS, bestMu_NLMS);

% Plot weights evolution
figure;
plot(w_hist_NLMS');
grid on;
title('NLMS Weights Evolution');
xlabel('Cycle (n)');
ylabel('Magnitude');
legend(arrayfun(@(i) sprintf('w_%d(n)', i - 1), 1:bestFilterOrder_NLMS, 'UniformOutput', false));

%% RLS: Plot time evolution of the filter weights for the best SNR configuration
[y_RLS, ~, w_hist_RLS] = rls_filter(noise, primary, bestFilterOrder_RLS, bestLambda_RLS);

% Plot weights evolution
figure;
plot(w_hist_RLS');
grid on;
title('RLS Weights Evolution');
xlabel('Cycle (n)');
ylabel('Magnitude');
legend(arrayfun(@(i) sprintf('w_%d(n)', i - 1), 1:bestFilterOrder_RLS, 'UniformOutput', false));