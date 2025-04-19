clear all; close all; clc;
addpath('C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Matlab\adaptivefilterMT\3_Baseline_Company_Data\Old unorganized stuff'); 


[audio, fs] = audioread("Separated_Signal_1.wav");

% Define parameters
%% windowLength is the size of each analysis window in samples.
windowLength = round(0.14*fs);  % window length is 140 ms 
%% The hop size (step size) determines how much the window moves forward.
hop_size = round(0.02*fs); % calculate next window by shifting 20ms worth of samples
%% overlap set the number of samples overlap between consecutive windows
overlap = windowLength-hop_size; 
%% % Number of Mel bands
numBands = 128; 
%% % Create a Hanning window
hannWin = hann(windowLength, 'periodic'); 

% Compute the Mel spectrogram
[s, f, t] = melSpectrogram(audio, fs, ...
                            "Window", kaiser(windowLength,18), ...
                           'OverlapLength', overlap, ...
                           'NumBands', numBands);

% Convert to dB scale for better visualization
s_db = 10 * log10(s + eps);

% Plot the Mel spectrogram
figure;
imagesc(t, f, s_db);
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Mel Spectrogram');
colorbar;
colormap magma; % Use 'jet' colormap for better visualization
set(gcf, 'Units', 'inches', 'Position', [1, 1, numel(audio)/fs, 1]);