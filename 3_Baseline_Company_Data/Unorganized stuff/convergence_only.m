clc;
clear;
close all;

%% Define the path for the data
[primary, fs] = audioread("C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Data from AI heathway\Data_ANC\Experiment_Data\Hospital Ambient Noises\NHS\1\primary.wav");   %noise + clean signal
[noise, ~] = audioread("C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Data from AI heathway\Data_ANC\Experiment_Data\Hospital Ambient Noises\NHS\1\secondary.wav");
[clean, ~] = audioread("C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Data from AI heathway\Data_ANC\Experiment_Data\Hospital Ambient Noises\NHS\1\ZCH0019.wav");


%% Define start time for trimming (1.7 seconds)
start_time = 1.7;  % in seconds
start_sample = round(start_time * fs);  % Convert time to sample index

%% Trim the signals to start from 1.7 seconds
primary = primary(start_sample:end);
noise = noise(start_sample:end);
clean = clean(start_sample:end);

% Find the minimum length
minLength = min([length(primary), length(noise), length(clean)]);

% Truncate signals to the same length
primary = primary(1:minLength);
noise = noise(1:minLength);
clean = clean(1:minLength);

%% Normalize signals to have similar range
primary = primary / (max(abs(primary)) + eps);
noise = noise / (max(abs(noise)) + eps);
clean = clean / (max(abs(clean)) + eps);

%% Bandpass Filtering to Remove Distortions
lowCutoff = 5; % Low cut frequency in Hz
highCutoff = 800; % High cut frequency in Hz
[b, a] = butter(4, [lowCutoff, highCutoff] / (fs/2), 'bandpass');
primary = filtfilt(b, a, primary);
noise = filtfilt(b, a, noise);

%% Check signal correlation
figure;
subplot(2,1,1);
[c1, lags] = xcorr(noise, primary-clean, 'normalized');
plot(lags, c1);
title('Correlation between reference noise and actual noise');
xlabel('Lag'); ylabel('Correlation');
grid on;

subplot(2,1,2);
[c2, lags] = xcorr(noise, clean, 'normalized');
plot(lags, c2);
title('Correlation between reference noise and clean signal');
xlabel('Lag'); ylabel('Correlation');
grid on;

% Calculate initial SNR - with error handling
noise_component = primary - clean;
if all(isfinite(clean)) && all(isfinite(noise_component)) && sum(noise_component.^2) > 0
    initial_SNR = snr(clean, noise_component);
    fprintf('Initial SNR: %.2f dB\n', initial_SNR);
else
    fprintf('Warning: Could not calculate initial SNR due to non-finite values or zero noise power\n');
    initial_SNR = NaN;
end

%% Parameters - ADJUSTED for better performance
mu_LMS = 0.05;      % Reduced from 0.66 to a more typical value
mu_NLMS = 0.3;      % Reduced from 0.96 to a more typical value
filterOrder = 32;   % Increased from 4 to better model real-world noise
lambda = 0.99;      % Corrected from 0.36 to a typical RLS forgetting factor
window_size = 5000; % Window size for MSE calculation

%% LMS Filter Implementation
% Initialize variables
error_LMS = zeros(1, length(primary));
output_LMS = zeros(1, length(primary));
w_LMS = zeros(filterOrder, 1);  % Filter weights initialization
MSE_LMS = zeros(1, length(primary));

% LMS filter loop
for n = filterOrder:length(primary)
    x = noise(n:-1:n-filterOrder+1);  % Input vector (reference noise)
    d = primary(n);  % Desired signal (primary = noise + clean)
    
    % Filter output and error calculation
    y_LMS = w_LMS' * x;        % Estimated noise
    e_LMS = d - y_LMS;         % Error = primary - estimated noise
    
    % Weight update
    w_LMS = w_LMS + mu_LMS * x * e_LMS;
    
    % Store results
    error_LMS(n) = e_LMS;
    output_LMS(n) = y_LMS;
    
    % Compute windowed MSE
    if n >= filterOrder + window_size
        MSE_LMS(n) = mean(error_LMS(n-window_size+1:n).^2);
    else
        MSE_LMS(n) = mean(error_LMS(filterOrder:n).^2);
    end
end

%% NLMS Filter Implementation
% Initialize variables
error_NLMS = zeros(1, length(primary));
output_NLMS = zeros(1, length(primary));
w_NLMS = zeros(filterOrder, 1);  % Filter weights initialization
MSE_NLMS = zeros(1, length(primary));

% NLMS filter loop
for n = filterOrder:length(primary)
    x = noise(n:-1:n-filterOrder+1);  % Input vector
    d = primary(n);  % Desired signal
    
    % Filter output and error calculation
    y_NLMS = w_NLMS' * x;      % Estimated noise
    e_NLMS = d - y_NLMS;       % Error = primary - estimated noise
    
    % Weight update with normalization
    norm_x = norm(x)^2 + 1e-6;  % Avoid division by zero
    w_NLMS = w_NLMS + (mu_NLMS / norm_x) * x * e_NLMS;
    
    % Store results
    error_NLMS(n) = e_NLMS;
    output_NLMS(n) = y_NLMS;
    
    % Compute windowed MSE
    if n >= filterOrder + window_size
        MSE_NLMS(n) = mean(error_NLMS(n-window_size+1:n).^2);
    else
        MSE_NLMS(n) = mean(error_NLMS(filterOrder:n).^2);
    end
end

%% RLS Filter Implementation
% Initialize RLS parameters
error_RLS = zeros(1, length(primary));
output_RLS = zeros(1, length(primary));
w_RLS = zeros(filterOrder, 1);  % Filter weights initialization
P = 1.0 * eye(filterOrder);  % Inverse correlation matrix - increased initial value
MSE_RLS = zeros(1, length(primary));
delta = 0.001;  % Regularization term

% RLS filter loop
for n = filterOrder:length(primary)
    x = noise(n:-1:n-filterOrder+1);  % Input vector (reference noise)
    d = primary(n);  % Desired signal (primary = noise + clean)
    
    % Filter output and error calculation
    y_RLS = w_RLS' * x;        % Estimated noise
    e_RLS = d - y_RLS;         % Error = primary - estimated noise
    
    % RLS update rule
    k = (P * x) / (lambda + x' * P * x);  % Kalman gain
    w_RLS = w_RLS + k * e_RLS;            % Weight update
    P = (P - k * x' * P) / lambda;        % Update inverse correlation matrix
    
    % Check for instability in P matrix and reset if necessary
    if any(isnan(P(:))) || any(isinf(P(:)))
        P = 100 * eye(filterOrder);
        fprintf('Warning: RLS P matrix reset at sample %d due to numerical instability\n', n);
    end
    
    % Store results
    error_RLS(n) = e_RLS;
    output_RLS(n) = y_RLS;
    
    % Compute windowed MSE
    if n >= filterOrder + window_size
        MSE_RLS(n) = mean(error_RLS(n-window_size+1:n).^2);
    else
        MSE_RLS(n) = mean(error_RLS(filterOrder:n).^2);
    end
end

%% Calculate final SNR for each method
% Recovered signals are primary minus estimated noise
recovered_LMS = primary - output_LMS';
recovered_NLMS = primary - output_NLMS';
recovered_RLS = primary - output_RLS';

% Replace any NaN or Inf values
recovered_LMS(~isfinite(recovered_LMS)) = 0;
recovered_NLMS(~isfinite(recovered_NLMS)) = 0;
recovered_RLS(~isfinite(recovered_RLS)) = 0;

% Calculate SNR with error handling
try
    error_signal_LMS = recovered_LMS - clean;
    error_signal_LMS(~isfinite(error_signal_LMS)) = 0;
    if sum(error_signal_LMS.^2) > 0
        SNR_LMS = snr(clean, error_signal_LMS);
        fprintf('Final SNR (LMS): %.2f dB\n', SNR_LMS);
    else
        SNR_LMS = Inf;
        fprintf('Warning: Cannot calculate LMS SNR (zero error power)\n');
    end
catch
    SNR_LMS = NaN;
    fprintf('Warning: Error calculating LMS SNR\n');
end

try
    error_signal_NLMS = recovered_NLMS - clean;
    error_signal_NLMS(~isfinite(error_signal_NLMS)) = 0;
    if sum(error_signal_NLMS.^2) > 0
        SNR_NLMS = snr(clean, error_signal_NLMS);
        fprintf('Final SNR (NLMS): %.2f dB\n', SNR_NLMS);
    else
        SNR_NLMS = Inf;
        fprintf('Warning: Cannot calculate NLMS SNR (zero error power)\n');
    end
catch
    SNR_NLMS = NaN;
    fprintf('Warning: Error calculating NLMS SNR\n');
end

try
    error_signal_RLS = recovered_RLS - clean;
    error_signal_RLS(~isfinite(error_signal_RLS)) = 0;
    if sum(error_signal_RLS.^2) > 0
        SNR_RLS = snr(clean, error_signal_RLS);
        fprintf('Final SNR (RLS): %.2f dB\n', SNR_RLS);
    else
        SNR_RLS = Inf;
        fprintf('Warning: Cannot calculate RLS SNR (zero error power)\n');
    end
catch
    SNR_RLS = NaN;
    fprintf('Warning: Error calculating RLS SNR\n');
end

%% Plot the Convergence
figure;

% Skip initial values which might have extreme peaks
skip_samples = filterOrder + 100;

% Clean up any Inf or NaN values for plotting
MSE_LMS(~isfinite(MSE_LMS)) = max(MSE_LMS(isfinite(MSE_LMS)));
MSE_NLMS(~isfinite(MSE_NLMS)) = max(MSE_NLMS(isfinite(MSE_NLMS)));
MSE_RLS(~isfinite(MSE_RLS)) = max(MSE_RLS(isfinite(MSE_RLS)));

% Subplot 1: LMS Filter Convergence
subplot(3,1,1);
semilogy(skip_samples:length(primary), MSE_LMS(skip_samples:end), 'r', 'LineWidth', 1.5);
grid on;
xlabel('Iteration (Samples)');
ylabel('MSE');
title('LMS Filter Convergence');
if sum(isfinite(MSE_LMS(skip_samples:end))) > 0
    ylim([min(MSE_LMS(isfinite(MSE_LMS(skip_samples:end))))*0.5, ...
          max(MSE_LMS(isfinite(MSE_LMS(skip_samples:end))))*2]);
end

% Subplot 2: NLMS Filter Convergence
subplot(3,1,2);
semilogy(skip_samples:length(primary), MSE_NLMS(skip_samples:end), 'b', 'LineWidth', 1.5);
grid on;
xlabel('Iteration (Samples)');
ylabel('MSE');
title('NLMS Filter Convergence');
if sum(isfinite(MSE_NLMS(skip_samples:end))) > 0
    ylim([min(MSE_NLMS(isfinite(MSE_NLMS(skip_samples:end))))*0.5, ...
          max(MSE_NLMS(isfinite(MSE_NLMS(skip_samples:end))))*2]);
end

% Subplot 3: RLS Filter Convergence
subplot(3,1,3);
semilogy(skip_samples:length(primary), MSE_RLS(skip_samples:end), 'g', 'LineWidth', 1.5);
grid on;
xlabel('Iteration (Samples)');
ylabel('MSE');
title('RLS Filter Convergence');
if sum(isfinite(MSE_RLS(skip_samples:end))) > 0
    ylim([min(MSE_RLS(isfinite(MSE_RLS(skip_samples:end))))*0.5, ...
          max(MSE_RLS(isfinite(MSE_RLS(skip_samples:end))))*2]);
end

% Add overall title
sgtitle('Convergence Comparison of Adaptive Filters');

%% Plot the time-domain signals for comparison
figure;

% Time vector
t = (0:length(primary)-1)/fs;

% Original signals
subplot(5,1,1);
plot(t, primary);
title('Primary Signal (Noisy)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(5,1,2);
plot(t, clean);
title('Clean Signal (Reference)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Recovered signals
subplot(5,1,3);
plot(t, recovered_LMS);
if isfinite(SNR_LMS)
    title(['LMS Recovered Signal (SNR: ', num2str(SNR_LMS, '%.2f'), ' dB)']);
else
    title('LMS Recovered Signal (SNR: N/A)');
end
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(5,1,4);
plot(t, recovered_NLMS);
if isfinite(SNR_NLMS)
    title(['NLMS Recovered Signal (SNR: ', num2str(SNR_NLMS, '%.2f'), ' dB)']);
else
    title('NLMS Recovered Signal (SNR: N/A)');
end
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(5,1,5);
plot(t, recovered_RLS);
if isfinite(SNR_RLS)
    title(['RLS Recovered Signal (SNR: ', num2str(SNR_RLS, '%.2f'), ' dB)']);
else
    title('RLS Recovered Signal (SNR: N/A)');
end
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

%% Optional: Plot filter weight evolution for each algorithm
% This can provide insight into how the filters are adapting
figure;

% Create subplots for first 4 weights of each filter
for i = 1:4
    subplot(3,4,(i-1)+1);
    plot((1:length(primary))/fs, [zeros(filterOrder-1,1); w_LMS(i)*ones(length(primary)-filterOrder+1,1)]);
    title(['LMS Weight ' num2str(i)]);
    xlabel('Time (s)');
    ylabel('Weight Value');
    grid on;
    
    subplot(3,4,(i-1)+5);
    plot((1:length(primary))/fs, [zeros(filterOrder-1,1); w_NLMS(i)*ones(length(primary)-filterOrder+1,1)]);
    title(['NLMS Weight ' num2str(i)]);
    xlabel('Time (s)');
    ylabel('Weight Value');
    grid on;
    
    subplot(3,4,(i-1)+9);
    plot((1:length(primary))/fs, [zeros(filterOrder-1,1); w_RLS(i)*ones(length(primary)-filterOrder+1,1)]);
    title(['RLS Weight ' num2str(i)]);
    xlabel('Time (s)');
    ylabel('Weight Value');
    grid on;
end

sgtitle('Evolution of Filter Weights');