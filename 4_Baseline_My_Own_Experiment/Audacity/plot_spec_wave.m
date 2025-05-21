% Load the WAV file
[inputSignal, fs] = audioread('C:\Users\eloma\Desktop\new_adaptivefilter\4_Baseline_My_Own_Experiment\Audacity\Clean Heartbeats\48kHz\Clean_tricuspid48khz_heartbeat_Jacob.wav');
inputSignal = inputSignal(:,1);  % Use first channel if stereo

% Name of audio type to have in figures and folder name
suffix = 'tricuspid48khz';

%Normalize  
inputSignal = inputSignal / max(abs(inputSignal));

% Parameters
windowLength = round(0.14 * fs);
hop_size = round(0.01 * fs);
overlap = windowLength - hop_size;
numBands = 256;

% Time vector for waveform
t_wave = ((0:length(inputSignal)-1) / fs);

% Compute Mel spectrogram
[s_mel, f, t_spec] = melSpectrogram(inputSignal, fs, ...
    "Window", kaiser(windowLength, 18), ...
    'OverlapLength', overlap, ...
    'NumBands', numBands);

s_mel_db = 10 * log10(s_mel + eps);
t_spec = t_spec - t_spec(1);  % Normalize time axis

% Plot
figure;
x_max = 10;

% 1. Plot waveform
ax1 = subplot(2,1,1);
plot(t_wave, inputSignal, 'k');
xlabel('Time (s)');
ylabel('Amplitude');
title('Waveform');
ylim([-1 1]);
xlim([0 x_max]);

% 2. Plot Mel spectrogram
ax2 = subplot(2,1,2);
imagesc(t_spec, f, s_mel_db);
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Mel Spectrogram');
colormap jet;
caxis([-60 max(s_mel_db(:))]);
ylim([0 2000]);
xlim([0 x_max]);

% Move colorbar to the right of the Mel spectrogram
cb = colorbar;
set(cb, 'Units', 'normalized');

% Get position of subplot (Mel spectrogram)
pos_ax2 = get(ax2, 'Position');
% Adjust colorbar to sit right next to the spectrogram
set(cb, 'Position', [pos_ax2(1) + pos_ax2(3) + 0.01, pos_ax2(2), 0.015, pos_ax2(4)]);

% Resize figure for saving
set(gcf, 'Units', 'inches', 'Position', [5.833333333333333,5.825,8.375,2.008333333333334]);
tightfig();
saveas(gcf, ['Wave_and_Mel_Spectrogram - ' suffix '.pdf']);

%% For comparing two sound files 

clear all;
close all;

% File paths
file1 = "C:\Users\eloma\Desktop\new_adaptivefilter\4_Baseline_My_Own_Experiment\Audacity\pre_clean_heart.wav";
file2 = "C:\Users\eloma\Desktop\new_adaptivefilter\4_Baseline_My_Own_Experiment\Audacity\post_clean_heart.wav";

% Suffixes for labeling and saving
suffix1 = 'Before processing';
suffix2 = 'After processing';

% Load audio
[input1, fs1] = audioread(file1);
[input2, fs2] = audioread(file2);

% Ensure mono and same sampling rate
input1 = input1(:,1);
input2 = input2(:,1);
if fs1 ~= fs2
    error('Sampling rates do not match!');
end
fs = fs1; % use common fs

% Normalize
input1 = input1 / max(abs(input1));
input2 = input2 / max(abs(input2));

% Mel Spectrogram Parameters
windowLength = round(0.14 * fs);
hop_size = round(0.01 * fs);
overlap = windowLength - hop_size;
numBands = 256;

% Time vectors
t_wave1 = (0:length(input1)-1) / fs;
t_wave2 = (0:length(input2)-1) / fs;

% Mel Spectrograms
[spec1, f, t1] = melSpectrogram(input1, fs, ...
    "Window", kaiser(windowLength, 18), ...
    "OverlapLength", overlap, ...
    "NumBands", numBands);
s1_db = 10 * log10(spec1 + eps);
t1 = t1 - t1(1);

[spec2, ~, t2] = melSpectrogram(input2, fs, ...
    "Window", kaiser(windowLength, 18), ...
    "OverlapLength", overlap, ...
    "NumBands", numBands);
s2_db = 10 * log10(spec2 + eps);
t2 = t2 - t2(1);

% Plotting
figure;

% Set x-limit
x_max = 10;

% 1. Waveform 1
ax1 = subplot(2,2,1);
plot(t_wave1, input1, 'k');
xlabel('Time (s)');
ylabel('Amplitude');
title(['Waveform - ' suffix1]);
ylim([-1 1]);
xlim([0 x_max]);

% 2. Waveform 2
ax2 = subplot(2,2,2);
plot(t_wave2, input2, 'k');
xlabel('Time (s)');
ylabel('Amplitude');
title(['Waveform - ' suffix2]);
ylim([-1 1]);
xlim([0 x_max]);

% 3. Mel Spectrogram 1
ax3 = subplot(2,2,3);
imagesc(t1, f, s1_db);
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title(['Mel Spectrogram - ' suffix1]);
colormap jet;
caxis([-60 max([s1_db(:); s2_db(:)])]);
ylim([0 2000]);
xlim([0 x_max]);

% 4. Mel Spectrogram 2
ax4 = subplot(2,2,4);
imagesc(t2, f, s2_db);
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title(['Mel Spectrogram - ' suffix2]);
colormap jet;
caxis([-60 max([s1_db(:); s2_db(:)])]);
ylim([0 2000]);
xlim([0 x_max]);

% Shared colorbar to the right of subplot(2,2,4)
cb = colorbar;
set(cb, 'Units', 'normalized');
pos_ax4 = get(ax4, 'Position');
set(cb, 'Position', [pos_ax4(1) + pos_ax4(3) + 0.01, pos_ax4(2), 0.015, pos_ax4(4)]);

% Resize figure
set(gcf, 'Units', 'inches', 'Position', [5.416666666666667,6.858333333333333,7.708333333333333,1.875000000000001]);
tightfig();

% Save figure
saveas(gcf, ['Wave_MelSpectrogram_' suffix1 '_vs_' suffix2 '.pdf']);

%% Helper function: tightfig
function tightfig()
    % Tighten the figure layout
    set(gcf, 'Units', 'Inches');
    pos = get(gcf, 'Position');
    set(gcf, 'PaperUnits', 'Inches');
    set(gcf, 'PaperSize', [pos(3) pos(4)]);
    set(gcf, 'PaperPosition', [0 0 pos(3) pos(4)]);
end


