close all
clear all

% Load the wav file
[heartbeat, Fs] = audioread('C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Data from AI heathway\Data_ANC\Clear_Normal_HS\ZCH0019.wav');
t = (0:length(heartbeat)-1) / Fs; % Time vector

% Load the S1 and S2 timestamps from the text file
data = readtable('C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Data from AI heathway\Data_ANC\Clear_Normal_HS\ZCH0019.txt', 'Format', '%f%f%s', 'Delimiter', '\t');

% Extract S1 and S2 timestamps
S1_times = data{strcmp(data.Var3, 'S1'), 1:2};
S2_times = data{strcmp(data.Var3, 'S2'), 1:2};

% Plot the heartbeat waveform
figure;
plot(t, heartbeat, 'b', 'HandleVisibility', 'off');
xlabel('Time (s)');
ylabel('Amplitude');
%xlim([1 5]);
ylim([-2 2]);
title('Heart Sound Waveform');
hold on;

% Define colors for S1 and S2
s1_color = [0 1 1];  % Cyan
s2_color = [1 0.5 0.5]; % Light Red

% Highlight S1 regions
for i = 1:size(S1_times, 1)
    x = [S1_times(i, 1) S1_times(i, 2)];
    y = [min(heartbeat) min(heartbeat) max(heartbeat) max(heartbeat)];
    patch(x([1 2 2 1]), y, s1_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

% Highlight S2 regions
for i = 1:size(S2_times, 1)
    x = [S2_times(i, 1) S2_times(i, 2)];
    y = [min(heartbeat) min(heartbeat) max(heartbeat) max(heartbeat)];
    patch(x([1 2 2 1]), y, s2_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

% Dummy plots for legend
plot(NaN, NaN, 'Color', s1_color, 'LineWidth', 5, 'DisplayName', 'S1');
plot(NaN, NaN, 'Color', s2_color, 'LineWidth', 5, 'DisplayName', 'S2');

legend;
hold off;


%%  Cropped example of normal heartbeat
close all
clear all

% Load the MP3 file
[heartbeat, Fs] = audioread('C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Data from AI heathway\Data_ANC\Clear_Normal_HS\ZCH0019.wav');
t = (0:length(heartbeat)-1) / Fs; % Time vector

% Define the starting time (10 seconds)
start_time = 10.1;

% Find the index corresponding to the start time
start_index = find(t >= start_time, 1);

% Trim the time and heartbeat signals from the start time
t_trimmed = t(start_index:end) - start_time; % Reset time after 10 sec
heartbeat_trimmed = heartbeat(start_index:end);

% Load the S1 and S2 timestamps from the text file
data = readtable('C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Data from AI heathway\Data_ANC\Clear_Normal_HS\ZCH0019.txt', 'Format', '%f%f%s', 'Delimiter', '\t');

% Extract S1 and S2 timestamps
S1_times = data{strcmp(data.Var3, 'S1'), 1:2};
S2_times = data{strcmp(data.Var3, 'S2'), 1:2};

% Create a figure with subplots
figure('Position', [100, 100, 1000, 500]);

% Subplot 1: Heartbeat Waveform
subplot(2,1,1);
plot(t_trimmed, heartbeat_trimmed, 'b', 'HandleVisibility', 'off');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0 2.5]);
ylim([-2 2]);
%title('Heartbeat Waveform');
hold on;

% Define colors for S1 and S2
s1_color = [0 1 1];  % Cyan
s2_color = [1 0.5 0.5]; % Light Red

% Highlight S1 regions
for i = 1:size(S1_times, 1)
    % Adjust the S1 times to the trimmed time vector
    S1_start = S1_times(i, 1) - start_time;
    S1_end = S1_times(i, 2) - start_time;
    x = [S1_start S1_end];
    y = [min(heartbeat_trimmed) min(heartbeat_trimmed) max(heartbeat_trimmed) max(heartbeat_trimmed)];
    patch(x([1 2 2 1]), y, s1_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

% Highlight S2 regions
for i = 1:size(S2_times, 1)
    % Adjust the S2 times to the trimmed time vector
    S2_start = S2_times(i, 1) - start_time;
    S2_end = S2_times(i, 2) - start_time;
    x = [S2_start S2_end];
    y = [min(heartbeat_trimmed) min(heartbeat_trimmed) max(heartbeat_trimmed) max(heartbeat_trimmed)];
    patch(x([1 2 2 1]), y, s2_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

% Dummy plots for legend
plot(NaN, NaN, 'Color', s1_color, 'LineWidth', 5, 'DisplayName', 'S1');
plot(NaN, NaN, 'Color', s2_color, 'LineWidth', 5, 'DisplayName', 'S2');

legend;
hold off;

% ---- Plot the Spectrogram (Linear Frequency) ----
subplot(2,1,2);
window = hann(256);  % Hann window of size
overlap = 200;       % Overlap length
nfft = 1024;         % FFT points for higher frequency resolution

% Create the spectrogram
[S, F, T] = spectrogram(heartbeat_trimmed, window, overlap, nfft, Fs);

% Plot the magnitude of the spectrogram
dBmin = -20; % Set a minimum threshold (experiment with values like -50 to -80) (clips the lower-intensity noise, making the background truly black)
S_dB = 20*log10(abs(S)); % Convert to dB scale
S_dB(S_dB < dBmin) = dBmin; % Clip values below the threshold
imagesc(T, F, S_dB);

axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
%title('Spectrogram of Heartbeat Sounds');
xlim([0 2.5]);
ylim([0 1000]);

colormap hot; % hot or turbo
colorbar;

% Crop the figure and save as PDF
tightfig();
saveas(gcf, 'heartbeat_plot.pdf');
%%

close all;
clear all;

% Load the lung sound WAV file
[lungs, Fs] = audioread('C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Data from AI heathway\Data_ANC\Clear_Normal_LS\steth_20190622_14_47_24_ok.wav');
t = (0:length(lungs)-1) / Fs; % Time vector

% Load the I and E timestamps from the text file
data = readtable('C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Data from AI heathway\Data_ANC\Clear_Normal_LS\steth_20190622_14_47_24.txt', 'Format', '%f%f%s', 'Delimiter', '\t');

% Extract I and E timestamps
I_times = data{strcmp(data.Var3, 'I'), 1:2};
E_times = data{strcmp(data.Var3, 'E'), 1:2};

% Create a figure with subplots
figure('Position', [100, 100, 1000, 500]);

% Subplot 1: Lung Sound Waveform
subplot(2,1,1);
plot(t, lungs, 'b', 'HandleVisibility', 'off');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0 10]);
ylim([-0.5 0.5]);
%title('Lung Sound Waveform');
hold on;

% Define colors for I and E
i_color = [0 1 0];  % Green for Inhalation (I)
e_color = [1 0 0];  % Red for Exhalation (E)

% Highlight I (Inhalation) regions
for i = 1:size(I_times, 1)
    x = [I_times(i, 1) I_times(i, 2)];
    y = [min(lungs) min(lungs) max(lungs) max(lungs)];
    patch(x([1 2 2 1]), y, i_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

% Highlight E (Exhalation) regions
for i = 1:size(E_times, 1)
    x = [E_times(i, 1) E_times(i, 2)];
    y = [min(lungs) min(lungs) max(lungs) max(lungs)];
    patch(x([1 2 2 1]), y, e_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

% Dummy plots for legend
plot(NaN, NaN, 'Color', i_color, 'LineWidth', 5, 'DisplayName', 'Inhalation (I)');
plot(NaN, NaN, 'Color', e_color, 'LineWidth', 5, 'DisplayName', 'Exhalation (E)');

legend;
hold off;

% ---- Plot the Spectrogram (Linear Frequency) ----
subplot(2,1,2);
window = hann(256);  % Hann window of size
overlap = 200;       % Overlap length
nfft = 1024;         % FFT points for higher frequency resolution

% Create the spectrogram
[S, F, T] = spectrogram(lungs, window, overlap, nfft, Fs);

% Plot the magnitude of the spectrogram
dBmin = -50; % Set a minimum threshold (experiment with values like -70 or -80)
S_dB = 20*log10(abs(S)); % Convert to dB scale
S_dB(S_dB < dBmin) = dBmin; % Clip values below the threshold
imagesc(T, F, S_dB);

axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
%title('Spectrogram of Lung Sounds');
xlim([0 10]);
ylim([0 4000]);

% % Hold on to add lines on top of the spectrogram
% hold on;
% 
% % Plot vertical lines to mark I (Inhalation) regions on the spectrogram
% for i = 1:size(I_times, 1)
%     % Convert time to indices for the spectrogram
%     [~, idx_start] = min(abs(T - I_times(i, 1)));
%     [~, idx_end] = min(abs(T - I_times(i, 2)));
%     % Plot vertical lines at the I boundaries
%     plot([T(idx_start) T(idx_start)], [2000 3000], 'Color', i_color, 'LineWidth', 2, 'HandleVisibility', 'off');
%     plot([T(idx_end) T(idx_end)], [2000 3000], 'Color', i_color, 'LineWidth', 2, 'HandleVisibility', 'off');
% end
% 
% % Plot vertical lines to mark E (Exhalation) regions on the spectrogram
% for i = 1:size(E_times, 1)
%     % Convert time to indices for the spectrogram
%     [~, idx_start] = min(abs(T - E_times(i, 1)));
%     [~, idx_end] = min(abs(T - E_times(i, 2)));
%     % Plot vertical lines at the E boundaries
%     plot([T(idx_start) T(idx_start)], [2000 3000], 'Color', e_color, 'LineWidth', 2, 'HandleVisibility', 'off');
%     plot([T(idx_end) T(idx_end)], [2000 3000], 'Color', e_color, 'LineWidth', 2, 'HandleVisibility', 'off');
% end

colormap hot; % hot or turbo
colorbar;

% Crop the figure and save as PDF
tightfig();
saveas(gcf, 'lung_plot.pdf');

function tightfig()
    % Tighten the figure by removing excess whitespace
    set(gcf, 'Units', 'Inches');
    pos = get(gcf, 'Position');
    set(gcf, 'PaperUnits', 'Inches');
    set(gcf, 'PaperSize', [pos(3) pos(4)]);
    set(gcf, 'PaperPosition', [0 0 pos(3) pos(4)]);
end