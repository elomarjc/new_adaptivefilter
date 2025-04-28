%% Clear workspace and close figures
clear all;
close all;
clc;

%% Load the primary and secondary signals (noise and reference)
[u, fs] = audioread("C:\Users\eloma\Desktop\new_adaptivefilter\3_Baseline_Company_Data\ICA\Experiment_Data_vs2 - ICA\Hospital Ambient Noises\NHS\1\primary.wav");
[d, ~] = audioread("C:\Users\eloma\Desktop\new_adaptivefilter\3_Baseline_Company_Data\ICA\Experiment_Data_vs2 - ICA\Hospital Ambient Noises\NHS\1\secondary.wav");
[x, ~] = audioread("C:\Users\eloma\Desktop\new_adaptivefilter\3_Baseline_Company_Data\ICA\Experiment_Data_vs2 - ICA\Hospital Ambient Noises\NHS\1\ZCH0019.wav");

% Clean signal is not needed for ICA, its just for visual comparison. 

%% Define start time for trimming (x seconds)
start_time = 1.7;  % in seconds
start_sample = round(start_time * fs);  % Convert time to sample index

% Trim the signals to start from x seconds
u = u(start_sample:end);
d = d(start_sample:end);
x = x(start_sample:end);

% Find the minimum length
minLength = min([length(u), length(d), length(x)]);

%Normalize
x = x / max(abs(x));
u = u / max(abs(u));
d = d / max(abs(d));

%% Bandpass Filtering to Remove Distortions
lowCutoff = 5; % Low cut frequency in Hz
highCutoff = 800; % High cut frequency in Hz
[b, a] = butter(4, [lowCutoff, highCutoff] / (fs/2), 'bandpass');
u = filtfilt(b, a, u);
d = filtfilt(b, a, d);

%% Center the data (as required by PowerICA)
Y = [u'; d'];  % Create a 2-by-n mixture matrix
Y = bsxfun(@minus, Y, mean(Y, 2));  % Center the data

%% Whiten the data (PowerICA requires whitening)   (selecting the most adequate for the matrix, would be the idea)
[E, D] = eig(cov(Y', 1));  % Eigenvalue decomposition of the covariance
[Ds, ord] = sort(diag(D), 'descend');  % Sort by decreasing eigenvalue
E = E(:, ord(1:2));  % Take the first 2 eigenvectors (for the 2 sources)
lam = Ds(1:2);  % Take the first 2 eigenvalues
whiteningMatrix = diag(1 ./ sqrt(lam)) * E';  % Whitening matrix
dewhiteningMatrix = E * diag(sqrt(lam));  % Dewhitening matrix
X = whiteningMatrix * Y;  % Whiten the mixed signals

%% Apply PowerICA algorithm
W0 = orth(randn(2, 2));  % Initial random demixing matrix
nonlin = 'tanh';  % ICA nonlinearity
mode = 'serial';  % ICA computation mode
[W_est, flg] = PowerICA(X, nonlin, W0, mode);  % Perform PowerICA

%% Estimate the independent components (ICs)
S_est = W_est * X;  % Estimated independent components

%% Estimate the mixing matrix (up to sign and permutation ambiguities)
A_est = dewhiteningMatrix * W_est';
fprintf('The PowerICA estimate of A is\n');
disp(A_est);

%% Save data
% Use current date for folder and filename
date_str = datestr(now, 'yyyy-dd-mm');  % e.g., '2025-19-04'

suffix = 'hospital ambient';

% Construct folder and filename
foldername = ['PowerICA' date_str '_' suffix];
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

S1 = S1(:);  % Use as-is from ICA
S2 = S2(:);
  
%% Calculate SNR (optional) - Assuming you have a clean signal available for SNR calculation
% % Intial snr
% initial_SNR = 10 * log10(sum(x.^2) / sum((x - d).^2));

% % SNR calculation after ICA separation (use the clean signal for SNR)
% filtered_SNR_1 = 10 * log10(sum(x.^2) / sum((x - filtered_signal_1).^2));
% filtered_SNR_2 = 10 * log10(sum(x.^2) / sum((x - filtered_signal_2).^2));
% 
% fprintf('\n--- Results ---\n');
% 
% fprintf('Initial SNR: %.4f dB\n', initial_SNR);  % Print initial SNR
% 
% fprintf('S1 SNR: %.4f dB\n', filtered_SNR_1);
% fprintf('S1 SNR Improvment: %.4f dB\n\n', filtered_SNR_1-initial_SNR);
% 
% fprintf('S2 SNR: %.4f dB\n', filtered_SNR_2);
% fprintf('S2 SNR Improvment: %.4f dB\n', filtered_SNR_2-initial_SNR);

%% Save parameters 

% Define result text file path (same folder as saved .mat)
results_txt_path = fullfile(foldername, 'best_parameters_summary.txt');

% Open file for writing (overwrite mode)
fid = fopen(results_txt_path, 'w');

% Check if file opened successfully
if fid == -1
    error('Failed to create result text file.');
end

fprintf(fid, 'The PowerICA estimate of A is \n');
fprintf(fid, ' %.4f \n', A_est);

% % Summary Section
% fprintf(fid, '--- Absolute Best Parameters ---\n');
% 
% fprintf(fid, 'Initial SNR: %.4f dB\n', initial_SNR);  % Print initial SNR
% fprintf(fid, 'S1 SNR: %.4f dB\n', filtered_SNR_1);
% fprintf(fid, 'S1 SNR Improvment: %.4f dB\n\n', filtered_SNR_1-initial_SNR);
% fprintf(fid, 'S2 SNR: %.4f dB\n', filtered_SNR_2);
% fprintf(fid, 'S2 SNR Improvment: %.4f dB\n', filtered_SNR_2-initial_SNR);
% 
% fprintf(fid, '\n');
% 
% fprintf(fid, '---------------------------------\n');

% Close file
fclose(fid);

%% Plot Mel Spectrogram for LMS, NLMS, and RLS in one figure

windowLength = round(0.14*fs);  % 140 ms window length
hop_size = round(0.02*fs);      % 20 ms hop size
overlap = windowLength - hop_size; 
numBands = 128;                 % Number of Mel bands
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

% Compute Mel spectrogram for S1
[s_S1, ~, ~] = melSpectrogram(S1, fs, ...
                            "Window", kaiser(windowLength,18), ...
                            'OverlapLength', overlap, ...
                            'NumBands', numBands);

s_S1_db = 10 * log10(s_S1 + eps);

% Compute Mel spectrogram for S2
[s_S2, ~, ~] = melSpectrogram(S2, fs, ...
                            "Window", kaiser(windowLength,18), ...
                            'OverlapLength', overlap, ...
                            'NumBands', numBands);

s_S2_db = 10 * log10(s_S2 + eps);

% Plot Mel spectrograms in a single figure (5 subplots)
figure;

% Clean signal
subplot(4,1,1); 
imagesc(t, f, s_x_db); 
axis xy; 
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Clean Signal');
colorbar;
colormap jet;
ylim([0 4000]);
xlim([0 5]);

% Primary signal (heart+noise)
subplot(4,1,2);
imagesc(t, f, s_u_db); 
axis xy; 
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Primary Signal');
colorbar;
colormap jet;
ylim([0 4000]);
xlim([0 5]);

% Heart seperated
subplot(4,1,3); 
imagesc(t, f, s_S1_db); 
axis xy; 
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Separated Signal 1 (Heart sound)');
colorbar;
colormap jet;
ylim([0 4000]);
xlim([0 5]);

% Noise seperated
subplot(4,1,4);
imagesc(t, f, s_S2_db); 
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Separated Signal 1 (Noise)');
colorbar;
colormap jet;
ylim([0 4000]);
xlim([0 5]);

% Set figure size and position
set(gcf, 'Units', 'inches', 'Position', [0, 0, 8.266666666666666, 9.866666666666667]);


% Save the figure
tightfig();
saveas(gcf, fullfile(foldername, 'Mel Spectrogram - PowerICA - hopital ambient.pdf'));

%% Plot error signals

% Make sure signals are column vectors
% x = x(:);
% S1 = S1(:);
% S2 = S2(:);

figure;
subplot(4, 1, 1);
plot(x);
title('Clean Signal');
xlabel('Sample Number');
ylabel('Amplitude');
ylim([-1 1]);
xlim([0 minLength]);

subplot(4, 1, 2);
plot(u);
title('Noisy Signal');
xlabel('Sample Number');
ylabel('Amplitude');
ylim([-1 1]);
xlim([0 minLength]);

subplot(4, 1, 3); 
plot(x-S1);
title('Error signal - Separated Signal 1 (Heart sound)');
xlabel('Sample Number');
ylabel('Amplitude');
ylim([-1 1]);
xlim([0 minLength]);

subplot(4, 1, 4);
plot(x-S2);
title('Error signal -Separated Signal 2 (Noise)');
xlabel('Sample Number');
ylabel('Amplitude');
ylim([-1 1]);
xlim([0 minLength]);

% Save the figure
tightfig();
saveas(gcf, fullfile(foldername, 'PowerICA Separation Results - hopital ambient.pdf'));


% Calculate initial SNR
SNR_seperate_1 = 10 * log10(sum(x.^2) / sum((x - d).^2)); 
SNR_seperate_2 = 10 * log10(sum(x.^2) / sum((x - d).^2));  

%% Save Figure trimmed
function tightfig()
    set(gcf, 'Units', 'Inches');
    pos = get(gcf, 'Position');
    set(gcf, 'PaperUnits', 'Inches');
    set(gcf, 'PaperSize', [pos(3) pos(4)]);
    set(gcf, 'PaperPosition', [0 0 pos(3) pos(4)]);
end