%% Clear workspace and close figures
clear all;
close all;
clc;

%% Load the primary and secondary signals (noise and reference)
% [u, fs] = audioread("C:\Users\eloma\Desktop\new_adaptivefilter\3_Baseline_Company_Data\ICA\Experiment_Data_vs2 - ICA\Artifacts\NHS\2\primary.wav");
% [d, ~] = audioread("C:\Users\eloma\Desktop\new_adaptivefilter\3_Baseline_Company_Data\ICA\Experiment_Data_vs2 - ICA\Artifacts\NHS\2\secondary.wav");
% [x, ~] = audioread("C:\Users\eloma\Desktop\new_adaptivefilter\3_Baseline_Company_Data\ICA\Experiment_Data_vs2 - ICA\Artifacts\NHS\2\ZCH0066.wav");

[u, fs] = audioread("C:\Users\eloma\Desktop\new_adaptivefilter\4_Baseline_My_Own_Experiment\ICA - Matlab\NHS data - ICA - Jacob\ambient noise\Primary.wav");
[d, ~] = audioread("C:\Users\eloma\Desktop\new_adaptivefilter\4_Baseline_My_Own_Experiment\ICA - Matlab\NHS data - ICA - Jacob\ambient noise\Secondary.wav");
[x, ~] = audioread("C:\Users\eloma\Desktop\new_adaptivefilter\4_Baseline_My_Own_Experiment\ICA - Matlab\NHS data - ICA - Jacob\ambient noise\Clean.wav");

% Clean signal is not needed for ICA, its just for visual comparison. 

suffix = 'ambient noise';

% Define the duration of the segment to extract (x seconds)
segment_duration = 6;  % in seconds
segment_samples = segment_duration * fs;  % Convert duration to sample count

% Ensure that the recordings are long enough and extract the first x seconds from each
if length(u) >= segment_samples
    u = u(1:segment_samples);  % Extract first x seconds
else
    error('Primary signal is shorter than x seconds.');
end

if length(d) >= segment_samples
    d = d(1:segment_samples);  % Extract first x seconds
else
    error('Secondary signal is shorter than 6 seconds.');
end

if length(x) >= segment_samples
    x = x(1:segment_samples);  % Extract first x seconds
else
    error('Clean signal is shorter than x seconds.');
end

% Ensure all signals are the same size by trimming to the smallest length
min_len = min([length(u), length(d), length(x)]);
u = u(1:min_len);
d = d(1:min_len);
x = x(1:min_len);

%% Calculate initial SNR
initial_SNR = 10 * log10(sum(x.^2) / sum((x - d).^2));  % The same as 10 * log10(sum(x.^2) / sum((x - d).^2))

%Normalize
x = x / max(abs(x));
u = u / max(abs(u));
d = d / max(abs(d));

%%
% Run FastICA
[S_est, W, A] = fastICA(u, d, 2);

t = (0:length(u)-1)/fs;  % Assuming you have sampling frequency fs

%% Save data
% Use current date for folder and filename
date_str = datestr(now, 'yyyy-dd-mm');  % e.g., '2025-19-04'

% Construct folder and filename
foldername = ['FastICA_' date_str '_' suffix];
filename = fullfile(foldername, [foldername '.mat']);

% Create folder if it doesn't exist
if ~exist(foldername, 'dir')
    mkdir(foldername);
end

%% Export separated signals after ICA (optional)
% Normalize each separated signal to prevent clipping
S1 = S_est(1, :) / max(abs(S_est(1, :)));
S2 = S_est(2, :) / max(abs(S_est(2, :)));

audiowrite(fullfile(foldername, 'Separated_Signal_1.wav'), S1, fs);  % Separated first source (noise)
audiowrite(fullfile(foldername,'Separated_Signal_2.wav'), S2, fs);  % Separated second source (heart/lung sound)

fprintf('Separated signals saved as audio files.\n\n');

% Ensure signals are column vectors
x = x(:); u = u(:); S1 = S1(:); S2 = S2(:);

%% SNR

% Calculate SNR in dB after ICA for both separated signals
SNR_post_S1 = 10 * log10(sum(x.^2) / sum((x - S1).^2));
SNR_post_S2 = 10 * log10(sum(x.^2) / sum((x - S2).^2));

%% Display Initial and Post-SNR values
fprintf('Initial SNR (before ICA): %.2f dB\n', initial_SNR);
fprintf('Post SNR (after ICA) for Separated Signal 1: %.2f dB\n', SNR_post_S1);
fprintf('Post SNR (after ICA) for Separated Signal 2: %.2f dB\n', SNR_post_S2);

fprintf('SNR improvement Separated Signal 1: %.2f dB\n', SNR_post_S1-initial_SNR);
fprintf('SNR improvement Separated Signal 2: %.2f dB\n', SNR_post_S2-initial_SNR);

% Define result text file path (same folder as saved .mat)
results_txt_path = fullfile(foldername, 'SNR-results.txt');

% Open file for writing (overwrite mode)
fid = fopen(results_txt_path, 'w');

% Check if file opened successfully
if fid == -1
    error('Failed to create result text file.');
end

fprintf(fid, 'Initial SNR (before ICA): %.2f dB\n', initial_SNR);
fprintf(fid, 'Post SNR (after ICA) for Separated Signal 1: %.2f dB\n', SNR_post_S1);
fprintf(fid, 'Post SNR (after ICA) for Separated Signal 2: %.2f dB\n', SNR_post_S2);

fprintf(fid, 'SNR improvement Separated Signal 1: %.2f dB\n', SNR_post_S1-initial_SNR);
fprintf(fid, 'SNR improvement Separated Signal 2: %.2f dB\n', SNR_post_S2-initial_SNR);

% Close file
fclose(fid);

%% Mel Spectrogram Visualization (for FastICA)
windowLength = round(0.14 * fs);  % 140 ms
hop_size = round(0.01 * fs);      % 10 ms
overlap = windowLength - hop_size;
numBands = 256;

% Clean signal
[s_x, f, t] = melSpectrogram(x, fs, ...
                            "Window", kaiser(windowLength,18), ...
                            'OverlapLength', overlap, ...
                            'NumBands', numBands);
s_x_db = 10 * log10(s_x + eps);

% Noisy input
[s_u, ~, ~] = melSpectrogram(u, fs, ...
                            "Window", kaiser(windowLength,18), ...
                            'OverlapLength', overlap, ...
                            'NumBands', numBands);
s_u_db = 10 * log10(s_u + eps);

% Noisy input
[s_d, ~, ~] = melSpectrogram(d, fs, ...
                            "Window", kaiser(windowLength,18), ...
                            'OverlapLength', overlap, ...
                            'NumBands', numBands);
s_d_db = 10 * log10(s_d + eps);

% Separated signal 1
[s_S1, ~, ~] = melSpectrogram(S1, fs, ...
                            "Window", kaiser(windowLength,18), ...
                            'OverlapLength', overlap, ...
                            'NumBands', numBands);
s_S1_db = 10 * log10(s_S1 + eps);

% Separated signal 2
[s_S2, ~, ~] = melSpectrogram(S2, fs, ...
                            "Window", kaiser(windowLength,18), ...
                            'OverlapLength', overlap, ...
                            'NumBands', numBands);
s_S2_db = 10 * log10(s_S2 + eps);

% Plot
t = t - t(1);  % Reset time axis

% Find global min and max for color axis
cmin = -60;
cmax = max([max(s_x_db(:)), max(s_u_db(:)),max(s_d_db(:)), max(s_S1_db(:)), max(s_S2_db(:))]);

ymax = 2000;
xmax = 5;

figure;
subplot(5,1,1);
imagesc(t, f, s_x_db); 
axis xy;
xlabel('Time (s)'); 
ylabel('Frequency (Hz)'); 
title('Clean Signal'); 
colorbar; ylim([0 ymax]); 
caxis([cmin cmax]); % Set color axis
xlim([0 xmax]);

subplot(5,1,2);
imagesc(t, f, s_u_db); 
axis xy;
xlabel('Time (s)'); 
ylabel('Frequency (Hz)'); 
title('Noisy Signal 1 (Clean + Noise)'); 
colorbar; ylim([0 ymax]); 
caxis([cmin cmax]); % Set color axis
xlim([0 xmax]);

subplot(5,1,3);
imagesc(t, f, s_u_db); 
axis xy;
xlabel('Time (s)'); 
ylabel('Frequency (Hz)'); 
title('Noisy Signal 2 (Clean + Noise)'); 
colorbar; ylim([0 ymax]); 
caxis([cmin cmax]); % Set color axis
xlim([0 xmax]);

subplot(5,1,4);
imagesc(t, f, s_S1_db); 
axis xy;
xlabel('Time (s)'); 
ylabel('Frequency (Hz)'); 
title('Separated Signal 1'); 
colorbar; 
ylim([0 ymax]); 
caxis([cmin cmax]); % Set color axis
xlim([0 xmax]);

subplot(5,1,5);
imagesc(t, f, s_S2_db); 
axis xy;
xlabel('Time (s)'); 
ylabel('Frequency (Hz)'); 
title('Separated Signal 2'); 
colorbar; 
ylim([0 ymax]); 
caxis([cmin cmax]); % Set color axis
xlim([0 xmax]);

% Set figure size and export
set(gcf, 'Units', 'inches', 'Position', [0, 0, 8.27, 9.8]);
tightfig();
% Construct filename
saveas(gcf, fullfile(foldername, ['Mel Spectrogram - FastICA - ' suffix '.pdf']));

%% Error Signal Plotting
minLength = min([length(x), length(S1), length(S2)]);

figure;
subplot(5,1,1);
plot(x); title('Clean Signal'); xlabel('Sample'); ylabel('Amplitude');
xlim([0 minLength]); ylim([-1 1]);

subplot(5,1,2);
plot(u); title('Noisy Signal 1'); xlabel('Sample'); ylabel('Amplitude');
xlim([0 minLength]); ylim([-1 1]);

subplot(5,1,3);
plot(d); title('Noisy Signal 2'); xlabel('Sample'); ylabel('Amplitude');
xlim([0 minLength]); ylim([-1 1]);

subplot(5,1,4);
plot(x-S1); title('Error Signal (Separated Signal 1)'); xlabel('Sample'); ylabel('Amplitude');
xlim([0 minLength]); ylim([-1 1]);

subplot(5,1,5);
plot(x-S2); title('Error Signal (Separated Signal 2)'); xlabel('Sample'); ylabel('Amplitude');
xlim([0 minLength]); ylim([-1 1]);

tightfig();
saveas(gcf, fullfile(foldername, ['Error Signal - FastICA - ' suffix '.pdf']));

%% Functions

function [S_est, W, A_whitened] = fastICA(u, d, num_components)
% fastICA - Separate  signal from noise using FastICA.
% 
% Inputs:
%   u, d - Mixed signals (each a 1D vector, same length)
%   num_components - Number of independent components to extract (usually 2)
%
% Outputs:
%   S_est - Estimated source signals (each row is one component)
%   W - Estimated demixing matrix
%   A_whitened - Whitening matrix

% Stack mixed signals into matrix
X = [u(:)'; d(:)'];  % shape: [2 x N]
[m, N] = size(X);

% === Step 1: Centering ===
X_mean = mean(X, 2);
X_centered = X - X_mean;

% === Step 2: Whitening ===
covX = cov(X_centered');  % covariance matrix
[E, D] = eig(covX);        % E: eigenvectors, D: eigenvalues
D_inv_sqrt = diag(1 ./ sqrt(diag(D)));
A_whitened = D_inv_sqrt * E';  % Whitening matrix
X_white = A_whitened * X_centered;

% === Step 3: FastICA using logcosh nonlinearity ===
W = zeros(num_components, m);
max_iter = 1000;
tol = 1e-6;

for p = 1:num_components
    w = randn(m, 1);
    w = w / norm(w);
    
    for iter = 1:max_iter
        w_old = w;
        
        % g(u) = tanh(u), g'(u) = 1 - tanh(u)^2
        u_proj = w' * X_white;
        g_u = tanh(u_proj);
        g_u_prime = 1 - g_u .^ 2;
        
        w = (X_white * g_u') / N - mean(g_u_prime) * w;
        
        % Decorrelate from previous components (deflation)
        if p > 1
            w = w - W(1:p-1,:)' * (W(1:p-1,:) * w);
        end
        
        % Normalize
        w = w / norm(w);
        
        % Check for convergence
        if norm(w - w_old) < tol || norm(w + w_old) < tol
            break;
        end
    end
    
    W(p, :) = w';
end

% === Step 4: Estimate source signals ===
S_est = W * X_white;  % Each row = one independent component

end

function tightfig()
    set(gcf, 'Units', 'Inches');
    pos = get(gcf, 'Position');
    set(gcf, 'PaperUnits', 'Inches');
    set(gcf, 'PaperSize', [pos(3) pos(4)]);
    set(gcf, 'PaperPosition', [0 0 pos(3) pos(4)]);
end