%% Clear workspace and close figures
clear all;
close all;
clc;

%% Load the primary and secondary signals (noise and reference)
[u, fs] = audioread("C:\Users\eloma\Desktop\new_adaptivefilter\3_Baseline_Company_Data\ICA\Artifacts - ICA data with 10ms delay/mic_1.wav");
[d, ~] = audioread("C:\Users\eloma\Desktop\new_adaptivefilter\3_Baseline_Company_Data\ICA\Artifacts - ICA data with 10ms delay/mic_2.wav");

%[u, fs] = audioread("C:\Users\eloma\Desktop\new_adaptivefilter\4_Baseline_My_Own_Experiment\ICA - Matlab\Raw - ICA - Measurement\Speech\Primary.wav");
%[d, ~] = audioread("C:\Users\eloma\Desktop\new_adaptivefilter\4_Baseline_My_Own_Experiment\ICA - Matlab\Raw - ICA - Measurement\Speech\Secondary.wav");


% Clean signal is not needed for ICA, its just for visual comparison. 

suffix = 'Artifacts-delayed';

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

%Normalize
u = u / max(abs(u));
d = d / max(abs(d));

% Ensure all signals are the same size by trimming to the smallest length
min_len = min([length(u), length(d)]);
u = u(1:min_len);
d = d(1:min_len);

% % Remove phase/time delay between u and d using cross-correlation
% [xcorr_val, lags] = xcorr(u, d);
% [~, max_idx] = max(abs(xcorr_val));  % Find lag with maximum correlation
% time_delay = lags(max_idx);          % Best alignment lag
% 
% % Align d to u by correcting the delay
% if time_delay > 0
%     d_aligned = [zeros(time_delay, 1); d];              % delay d forward
%     d_aligned = d_aligned(1:length(d));                 % trim to original length
% elseif time_delay < 0
%     d_aligned = d(-time_delay + 1:end);                 % advance d
%     d_aligned = [d_aligned; zeros(-time_delay, 1)];     % pad to original length
% else
%     d_aligned = d;  % already aligned
% end
% 
% % Trim both signals to same length again after alignment
% min_len = min(length(u), length(d_aligned));
% u = u(1:min_len);
% d = d_aligned(1:min_len);

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
u = u(:); S1 = S1(:); S2 = S2(:);

%% Mel Spectrogram Visualization (for FastICA)
windowLength = round(0.14 * fs);  % 140 ms
hop_size = round(0.01 * fs);      % 10 ms
overlap = windowLength - hop_size;
numBands = 256;

% Noisy input
[s_u, f, t] = melSpectrogram(u, fs, ...
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
cmin = -50;
cmax = max([max(s_u_db(:)),max(s_d_db(:)), max(s_S1_db(:)), max(s_S2_db(:))]);

ymax = 2000;
xmax = 5;

figure;
subplot(4,1,1);
imagesc(t, f, s_u_db); 
axis xy;
xlabel('Time (s)'); 
ylabel('Frequency (Hz)'); 
title('Noisy Signal 1 (Clean + Noise)'); 
colorbar; ylim([0 ymax]); 
caxis([cmin cmax]); % Set color axis
xlim([0 xmax]);

subplot(4,1,2);
imagesc(t, f, s_u_db); 
axis xy;
xlabel('Time (s)'); 
ylabel('Frequency (Hz)'); 
title('Noisy Signal 2 (Clean + Noise)'); 
colorbar; ylim([0 ymax]); 
caxis([cmin cmax]); % Set color axis
xlim([0 xmax]);

subplot(4,1,3);
imagesc(t, f, s_S1_db); 
axis xy;
xlabel('Time (s)'); 
ylabel('Frequency (Hz)'); 
title('Separated Signal 1'); 
colorbar; 
ylim([0 ymax]); 
caxis([cmin cmax]); % Set color axis
xlim([0 xmax]);

subplot(4,1,4);
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

%% Waveform Plotting
figure;
subplot(4,1,1);
plot(u); title('Noisy Signal 1'); xlabel('Sample'); ylabel('Amplitude');
xlim([0 min_len]); ylim([-1 1]);

subplot(4,1,2);
plot(d); title('Noisy Signal 2'); xlabel('Sample'); ylabel('Amplitude');
xlim([0 min_len]); ylim([-1 1]);

subplot(4,1,3);
plot(S1); title('Separated Signal 1'); xlabel('Sample'); ylabel('Amplitude');
xlim([0 min_len]); ylim([-1 1]);

subplot(4,1,4);
plot(S2); title('Separated Signal 2'); xlabel('Sample'); ylabel('Amplitude');
xlim([0 min_len]); ylim([-1 1]);

tightfig();
saveas(gcf, fullfile(foldername, ['Waveform - FastICA - ' suffix '.pdf']));

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