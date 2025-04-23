clc;
clear;
close all;

%% Load real-world signals
[primary, fs]= audioread("C:\Users\eloma\Desktop\new_adaptivefilter\4_Baseline_My_Own_Experiment\NHS data\speech\Primary.wav");  
[reference, fs] = audioread("C:\Users\eloma\Desktop\new_adaptivefilter\4_Baseline_My_Own_Experiment\NHS data\speech\Secondary.wav");     

% Normalize
primary = primary / max(abs(primary));
reference = reference / max(abs(reference));

% Filter configs
filter_lengths = [64, 128, 256];
mu_vals = [0.01, 0.05, 0.1];          
lambda_vals = [0.98, 0.99];      

% Data storage
results = struct('method', {}, 'M', {}, 'param', {}, 'output', {}, 'weights', {}, 'est_snr', {}, 'error', {});

% Run LMS
for M = filter_lengths
    for mu = mu_vals
        w = zeros(M, 1);
        y = zeros(size(primary));
        e = zeros(size(primary));

        for n = M:length(primary)
            x = reference(n:-1:n-M+1);
            y(n) = w' * x;
            e(n) = primary(n) - y(n);
            w = w + mu * x * e(n);
        end

        est_snr = 10 * log10(var(y) / var(e)); 

        results(end+1) = struct('method', 'LMS', 'M', M, 'param', mu, ...
            'output', y, 'weights', w, 'est_snr', est_snr, 'error', e);
    end
end

% Run NLMS
for M = filter_lengths
    for mu = mu_vals
        w = zeros(M, 1);
        y = zeros(size(primary));
        e = zeros(size(primary));

        for n = M:length(primary)
            x = reference(n:-1:n-M+1);
            norm_factor = x' * x + 1e-6;
            y(n) = w' * x;
            e(n) = primary(n) - y(n);
            w = w + (mu / norm_factor) * x * e(n);
        end

        est_snr = 10 * log10(var(y) / var(e));

        results(end+1) = struct('method', 'NLMS', 'M', M, 'param', mu, ...
            'output', y, 'weights', w, 'est_snr', est_snr, 'error', e);
    end
end

% Run RLS
for M = filter_lengths
    for lambda = lambda_vals
        w = zeros(M, 1);
        P = eye(M) * 1000;
        y = zeros(size(primary));
        e = zeros(size(primary));

        for n = M:length(primary)
            x = reference(n:-1:n-M+1);
            pi_vec = P * x;
            k = pi_vec / (lambda + x' * pi_vec);
            y(n) = w' * x;
            e(n) = primary(n) - y(n);
            w = w + k * e(n);
            P = (P - k * x' * P) / lambda;
        end

        est_snr = 10 * log10(var(y) / var(e));

        results(end+1) = struct('method', 'RLS', 'M', M, 'param', lambda, ...
            'output', y, 'weights', w, 'est_snr', est_snr, 'error', e);
    end
end

%% Identify best configs per method
methods = {'LMS', 'NLMS', 'RLS'};
best_outputs = struct();

for m = 1:length(methods)
    method_results = results(strcmp({results.method}, methods{m}));
    [~, idx] = max([method_results.est_snr]);
    best_outputs.(methods{m}) = method_results(idx);
end

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
plot(best_outputs.LMS.output);
title('Error signal (LMS)');
xlabel('Sample Number');
ylabel('Amplitude');
ylim([-1.1 1.1]);
xlim([0 length(primary)]);

subplot(5, 1, 4);
plot(best_outputs.NLMS.output);
title('Error signal (NLMS)');
xlabel('Sample Number');
ylabel('Amplitude');
ylim([-1.1 1.1]);
xlim([0 length(primary)]);

subplot(5, 1, 5);
plot(best_outputs.RLS.output);
title('Error signal (RLS)');
xlabel('Sample Number');
ylabel('Amplitude');
ylim([-1.1 1.1]);
xlim([0 length(primary)]);


%% Plot frequency response of best filter
figure;
for i = 1:length(methods)
    subplot(length(methods), 1, i);
    w = best_outputs.(methods{i}).weights;
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
    'LMS', best_outputs.LMS.output;
    'NLMS', best_outputs.NLMS.output;
    'RLS', best_outputs.RLS.output;
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

%% Plot SNR as a function of filter order
figure;
hold on;
colors = lines(3);
markers = {'o', 's', '^'};
legend_entries = {};

for i = 1:length(methods)
    method = methods{i};
    method_results = results(strcmp({results.method}, method));

    % Group by filter order
    unique_orders = unique([method_results.M]);
    snr_per_order = zeros(size(unique_orders));
    
    for j = 1:length(unique_orders)
        order = unique_orders(j);
        order_results = method_results([method_results.M] == order);
        snr_per_order(j) = max([order_results.est_snr]);
    end

    plot(unique_orders, snr_per_order, ...
        ['-' markers{i}], ...
        'Color', colors(i,:), ...
        'LineWidth', 1.5, ...
        'MarkerSize', 8);

    legend_entries{end+1} = method;
end

title('Estimated SNR vs Filter Order');
xlabel('Filter Order (M)');
ylabel('Estimated SNR (dB)');
legend(legend_entries, 'Location', 'best');
grid on;
hold off;

