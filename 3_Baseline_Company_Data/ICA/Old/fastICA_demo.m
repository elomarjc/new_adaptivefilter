% fastICA_demo.m
% Demonstration script for FastICA heartbeat separation

[u, fs] = audioread("C:\Users\eloma\Desktop\new_adaptivefilter\3_Baseline_Company_Data\ICA\Experiment_Data_vs2 - ICA\Hospital Ambient Noises\NHS\1\primary.wav");
[d, ~] = audioread("C:\Users\eloma\Desktop\new_adaptivefilter\3_Baseline_Company_Data\ICA\Experiment_Data_vs2 - ICA\Hospital Ambient Noises\NHS\1\secondary.wav");


Y = [u'; d'];  % Create a 2-by-n mixture matrix

% Apply FastICA
[separated_heartbeat, separated_noise] = fastICA_heartbeat(Y, fs);

% Display correlation with original heartbeat
corr_val = abs(corr(heartbeat, separated_heartbeat));
fprintf('Correlation with original heartbeat: %.4f\n', corr_val);

% Compare original and extracted signals
figure;
subplot(2,1,1);
plot(t, heartbeat);
title('Original Heartbeat Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(t, separated_heartbeat);
title('Extracted Heartbeat Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;


%%

% Demonstrates the use of fastICA_heartbeat with simulated data
    
% Parameters
fs = 1000;              % Sampling frequency (Hz)
duration = 10;          % Duration (seconds)
t = (0:1/fs:duration-1/fs)';  % Time vector

% Generate a simulated heartbeat signal (roughly 60 BPM)
heart_rate = 60/60;     % Heart rate in Hz (60 BPM)
heartbeat = sin(2*pi*heart_rate*t) .* exp(-((mod(t, 1/heart_rate)-0.2).^2)/0.01);
heartbeat = heartbeat / std(heartbeat);

% Add noise
noise1 = 0.5*sin(2*pi*50*t);  % Power line interference
noise2 = 0.3*randn(size(t));  % Random noise

% Create two mixed signals
mixing_matrix = [0.8 0.2; 0.3 0.7];
mixed_signal = [heartbeat + noise1, heartbeat + noise2] * mixing_matrix;

% Apply FastICA
[separated_heartbeat, separated_noise] = fastICA_heartbeat(mixed_signal, fs);

% Display correlation with original heartbeat
corr_val = abs(corr(heartbeat, separated_heartbeat));
fprintf('Correlation with original heartbeat: %.4f\n', corr_val);

% Compare original and extracted signals
figure;
subplot(2,1,1);
plot(t, heartbeat);
title('Original Heartbeat Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(t, separated_heartbeat);
title('Extracted Heartbeat Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;