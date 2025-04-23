clc;
clear;
close all;

%% Load real-world signals
[primary, fs] = audioread("C:\Users\eloma\Desktop\new_adaptivefilter\4_Baseline_My_Own_Experiment\NHS data\speech\Primary.wav");  
[reference, fs] = audioread("C:\Users\eloma\Desktop\new_adaptivefilter\4_Baseline_My_Own_Experiment\NHS data\speech\Secondary.wav");     

% Normalize
primary = primary / max(abs(primary));
reference = reference / max(abs(reference));

%% Set parameters manually
M = 128;              % Filter length
mu_LMS = 0.01;        % LMS step size
mu_NLMS = 0.05;       % NLMS step size
lambda_RLS = 0.99;    % RLS forgetting factor

%% Data storage
results = struct();

%% Run LMS
w = zeros(M, 1);
y = zeros(size(primary));
e = zeros(size(primary));

for n = M:length(primary)
    x = reference(n:-1:n-M+1);
    y(n) = w' * x;
    e(n) = primary(n) - y(n);
    w = w + mu_LMS * x * e(n);
end

est_snr = 10 * log10(var(y) / var(e));

results.LMS = struct('output', y, 'weights', w, 'est_snr', est_snr, 'error', e);

%% Run NLMS
w = zeros(M, 1);
y = zeros(size(primary));
e = zeros(size(primary));

for n = M:length(primary)
    x = reference(n:-1:n-M+1);
    norm_factor = x' * x + 1e-6;
    y(n) = w' * x;
    e(n) = primary(n) - y(n);
    w = w + (mu_NLMS / norm_factor) * x * e(n);
end

est_snr = 10 * log10(var(y) / var(e));

results.NLMS = struct('output', y, 'weights', w, 'est_snr', est_snr, 'error', e);

%% Run RLS
w = zeros(M, 1);
P = eye(M) * 1000;
y = zeros(size(primary));
e = zeros(size(primary));

for n = M:length(primary)
    x = reference(n:-1:n-M+1);
    pi_vec = P * x;
    k = pi_vec / (lambda_RLS + x' * pi_vec);
    y(n) = w' * x;
    e(n) = primary(n) - y(n);
    w = w + k * e(n);
    P = (P - k * x' * P) / lambda_RLS;
end

est_snr = 10 * log10(var(y) / var(e));

results.RLS = struct('output', y, 'weights', w, 'est_snr', est_snr, 'error', e);

%% Plot error signals
figure;
subplot(5, 1, 1);
plot(primary);
title('Clean Signal');
xlabel('Sample Number');
ylabel('Amplitude');
ylim([-1.1 1.1]);
xlim([0 length(primary)]);

subplot(5, 1, 2);
plot(reference);
title('Noisy Signal');
xlabel('Sample Number');
ylabel('Amplitude');
ylim([-1.1 1.1]);
xlim([0 length(primary)]);

subplot(5, 1, 3); 
plot(results.LMS.output);
title('Error signal (LMS)');
xlabel('Sample Number');
ylabel('Amplitude');
ylim([-1.1 1.1]);
xlim([0 length(primary)]);

subplot(5, 1, 4);
plot(results.NLMS.output);
title('Error signal (NLMS)');
xlabel('Sample Number');
ylabel('Amplitude');
ylim([-1.1 1.1]);
xlim([0 length(primary)]);

subplot(5, 1, 5);
plot(results.RLS.output);
title('Error signal (RLS)');
xlabel('Sample Number');
ylabel('Amplitude');
ylim([-1.1 1.1]);
xlim([0 length(primary)]);

%% Plot frequency response of best filter
methods = {'LMS', 'NLMS', 'RLS'};

figure;
for i = 1:length(methods)
    subplot(length(methods), 1, i);
    w = results.(methods{i}).weights;
    [H, F] = freqz(w, 1, 1024, fs);
    plot(F, 20*log10(abs(H)), 'LineWidth', 1.5);
    title(sprintf('Frequency Response (%s)', methods{i}));
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    grid on;
end

%% Mel Spectrograms
windowLength = round(0.14 * fs);
hop_size = round(0.02 * fs);
overlap = windowLength - hop_size;
numBands = 128;

hannWin = kaiser(windowLength, 18);

signals_to_plot = {
    'Primary', primary;
    'Reference', reference;
    'LMS', results.LMS.output;
    'NLMS', results.NLMS.output;
    'RLS', results.RLS.output;
};

figure;
for i = 1:size(signals_to_plot,1)
    subplot(5,1,i);
    [S, f, t] = melSpectrogram(signals_to_plot{i,2}, fs, ...
        'Window', hannWin, ...
        'OverlapLength', overlap, ...
        'NumBands', numBands);
    imagesc(t, f, 10*log10(S + eps));
    axis xy;
    title([signals_to_plot{i,1} ' Signal']);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    ylim([0 4000]);
    colorbar;
end
colormap jet;

%% Save audio
output_folder = 'Filtered_Results_newthing_baseline';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

audiowrite(fullfile(output_folder, 'Clean_Primary.wav'), primary, fs);
audiowrite(fullfile(output_folder, 'Noisy_Reference.wav'), reference, fs);
audiowrite(fullfile(output_folder, 'Filtered_LMS.wav'), results.LMS.output, fs);
audiowrite(fullfile(output_folder, 'Filtered_NLMS.wav'), results.NLMS.output, fs);
audiowrite(fullfile(output_folder, 'Filtered_RLS.wav'), results.RLS.output, fs);
