clc;
clear;
close all;

%% Define the path for the data
[u, fs] = audioread("C:\Users\eloma\Desktop\new_adaptivefilter\4_Baseline_My_Own_Experiment\NHS data\speech\Primary.wav");  % Primary - Noisy signal with heartbeat   
[d, ~] = audioread("C:\Users\eloma\Desktop\new_adaptivefilter\4_Baseline_My_Own_Experiment\NHS data\speech\Secondary.wav"); % Reference noise signal     


% Name of audio type to have in figures and folder name
suffix = 'speech';

%% Parameters (Single values)
mu_LMS = 0.3;       % For LMS
mu_NLMS = 0.03;      % For NLMS
lambda_RLS = 0.9999;  % For RLS

M_LMS = 8;         % LMS filter length
M_NLMS = 8;        % NLMS filter length
M_RLS = 10;        % RLS filter length

%% Filtering (using parfor)
filteredSignals = cell(3,1);

parfor i = 1:3
    switch i
        case 1 % LMS
            [y, ~] = lms_filter(u, d, M_LMS, mu_LMS);
            filteredSignals{i} = u - y;
        case 2 % NLMS
            [y, ~] = nlms_filter(u, d, M_NLMS, mu_NLMS);
            filteredSignals{i} = u - y;
        case 3 % RLS
            [y, ~] = rls_filter(u, d, M_RLS, lambda_RLS);
            filteredSignals{i} = u - y;
    end
end

filtered_LMS = filteredSignals{1};
filtered_NLMS = filteredSignals{2};
filtered_RLS = filteredSignals{3};

%% Normalize signals for saving
filtered_LMS = filtered_LMS / max(abs(filtered_LMS));
filtered_NLMS = filtered_NLMS / max(abs(filtered_NLMS));
filtered_RLS = filtered_RLS / max(abs(filtered_RLS));

%% Save Audio
date_str = datestr(now, 'yyyy-dd-mm');  % e.g., '2025-19-04'
foldername = [date_str '_' suffix]; % Construct folder

% Create folder if it doesn't exist
if ~exist(foldername, 'dir')
    mkdir(foldername);
end

audiowrite(fullfile(foldername, 'Filtered_LMS.wav'), filtered_LMS, fs);
audiowrite(fullfile(foldername, 'Filtered_NLMS.wav'), filtered_NLMS, fs);
audiowrite(fullfile(foldername, 'Filtered_RLS.wav'), filtered_RLS, fs);

fprintf('Filtered signals saved successfully!\n');

%% Plot Mel Spectrogram for LMS, NLMS, and RLS in one figure

windowLength = round(0.14*fs);  % 140 ms window length
hop_size = round(0.01*fs);      % 10 ms hop size
overlap = windowLength - hop_size; 
numBands = 256;                 % Number of Mel bands
hannWin = hann(windowLength, 'periodic'); 


% Compute Mel spectrogram for primary signal (heart+noise)
[s_u, f, t] = melSpectrogram(u, fs, ...
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
cmax = max([max(s_u_db(:)), max(s_filtered_LMS_db(:)), max(s_filtered_NLMS_db(:)), max(s_filtered_RLS_db(:))]);

% Plot Mel spectrograms in a single figure (4 subplots)
figure;

% Primary signal (heart+noise)
subplot(4,1,1);
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
subplot(4,1,2); 
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
subplot(4,1,3);
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
subplot(4,1,4);
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

% %% Plot Filtered signals
% 
% figure;
% subplot(4, 1, 1);
% plot(u);
% title('Primary Signal (Clean + Noise)');
% xlabel('Sample Number');
% ylabel('Amplitude');
% ylim([-1.1 1.1]);
% xlim([0 length(u)]);
% 
% subplot(4, 1, 2); 
% plot(filtered_LMS);
% title('Filtered signal (LMS)');
% xlabel('Sample Number');
% ylabel('Amplitude');
% ylim([-1.1 1.1]);
% xlim([0 length(u)]);
% 
% subplot(4, 1, 3);
% plot(filtered_NLMS);
% title('Filtered signal (NLMS)');
% xlabel('Sample Number');
% ylabel('Amplitude');
% ylim([-1.1 1.1]);
% xlim([0 length(u)]);
% 
% subplot(4, 1, 4);
% plot(filtered_RLS);
% title('Filtered signal (RLS)');
% xlabel('Sample Number');
% ylabel('Amplitude');
% ylim([-1.1 1.1]);
% xlim([0 length(u)]);
% 
% % Save figure
% tightfig();
% saveas(gcf, fullfile(foldername, ['Filtered signals - ' suffix '.pdf']));

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

        %w_hist(:, n) = w;
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

        %w_hist(:, n) = w;
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

        %w_hist(:, n) = w;
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