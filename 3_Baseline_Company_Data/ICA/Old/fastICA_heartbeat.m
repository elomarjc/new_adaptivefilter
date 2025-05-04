function [heartbeat_signal, noise] = fastICA_heartbeat(mixed_signal, fs, visualize)
% FASTICA_HEARTBEAT - Separate heartbeat signal from noise using FastICA
%
% Inputs:
%   mixed_signal - Matrix where each column is a different recording channel
%   fs - Sampling frequency in Hz
%   visualize - Boolean flag to enable/disable visualization (default: true)
%
% Outputs:
%   heartbeat_signal - Extracted heartbeat signal
%   noise - Extracted noise component
%
% Implementation based on the FastICA algorithm

% Set default for visualization parameter
if nargin < 3
    visualize = true;
end

%% Step 1: Preprocessing
% Center the data (subtract mean)
X = mixed_signal - mean(mixed_signal, 1);

% Step 2: Whitening
% Calculate covariance matrix
C = X' * X / (size(X, 1) - 1);

% Eigenvalue decomposition
[E, D] = eig(C);
D = diag(D);

% Sort eigenvalues in descending order
[D, order] = sort(D, 'descend');
E = E(:, order);

% Whitening transformation
X_whitened = X * E * diag(1./sqrt(D));

%% Step 3: FastICA
% Number of components
n_comp = size(X, 2);

% Initialize mixing matrix
W = randn(n_comp, n_comp);

% Orthogonalize initial weight matrix
W = orth(W);

% Maximum number of iterations
max_iter = 1000;
tolerance = 1e-6;

% FastICA with symmetric orthogonalization
for iter = 1:max_iter
    W_old = W;
    
    % Calculate G and G'
    % Using G(y) = log(cosh(y)) as a robust non-quadratic function
    Y = X_whitened * W;
    
    % g(y) = tanh(y) - derivative of log(cosh(y))
    gY = tanh(Y);
    
    % g'(y) = 1 - tanh^2(y)
    gPrimeY = 1 - tanh(Y).^2;
    
    % Update W
    W = (X_whitened' * gY) / size(X_whitened, 1) - ...
        mean(gPrimeY, 1)' * W;
    
    % Symmetric orthogonalization
    W = W * (W' * W)^(-0.5);
    
    % Check for convergence
    if norm(abs(W) - abs(W_old)) < tolerance
        break;
    end
end

% Extract independent components
S = X_whitened * W;

%% Step 4: Component selection (identify the heartbeat signal)
% For heartbeat signals, we can use several criteria:
% 1. Periodicity analysis using autocorrelation
% 2. Frequency domain analysis
% We'll use a simple frequency domain approach

% Calculate power spectrum of each component
nfft = 2^nextpow2(size(S, 1));
all_psd = zeros(floor(nfft/2)+1, n_comp);

for i = 1:n_comp
    [pxx, f] = periodogram(S(:,i), hamming(length(S(:,i))), nfft, fs);
    all_psd(:, i) = pxx;
end

% Find frequency range for typical heart rates (0.5-3 Hz, or 30-180 BPM)
freq_range = [0.5 3.0];
freq_idx = find(f >= freq_range(1) & f <= freq_range(2));

% Calculate power in heartbeat frequency range for each component
power_in_hr_range = zeros(n_comp, 1);
for i = 1:n_comp
    power_in_hr_range(i) = sum(all_psd(freq_idx, i));
end

% Component with highest power in heartbeat frequency range
[~, heartbeat_idx] = max(power_in_hr_range);

% Get heartbeat and noise components
heartbeat_signal = S(:, heartbeat_idx);
noise_idx = setdiff(1:n_comp, heartbeat_idx);
noise = S(:, noise_idx);

% Normalize the heartbeat signal
heartbeat_signal = (heartbeat_signal - mean(heartbeat_signal)) / std(heartbeat_signal);

%% Visualization
if visualize
    % Time domain plot
    t = (0:size(mixed_signal,1)-1)' / fs;
    
    figure('Name', 'FastICA for Heartbeat Separation');
    
    % Original mixed signals
    subplot(3, 1, 1);
    plot(t, X);
    title('Original Mixed Signals');
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
    
    % Extracted heartbeat
    subplot(3, 1, 2);
    plot(t, heartbeat_signal);
    title('Extracted Heartbeat Signal');
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
    
    % Frequency domain analysis
    subplot(3, 1, 3);
    [pxx, f] = periodogram(heartbeat_signal, hamming(length(heartbeat_signal)), nfft, fs);
    plot(f, 10*log10(pxx));
    title('Power Spectrum of Heartbeat Signal');
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    xlim([0 5]); % Focusing on 0-5 Hz where heart rate typically falls
    grid on;
    
    % Mark heart rate frequency range
    hold on;
    y_lim = get(gca, 'YLim');
    area([freq_range(1) freq_range(2)], [y_lim(2) y_lim(2)], 'BaseValue', y_lim(1), ...
         'FaceColor', [0.8 0.8 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    text(mean(freq_range), y_lim(1) + 0.9*(y_lim(2)-y_lim(1)), 'Heart Rate Range', ...
         'HorizontalAlignment', 'center');
    hold off;
end
end

% Example usage:
function demo_fastICA_heartbeat()
    % This function demonstrates the use of fastICA_heartbeat
    % with simulated data
    
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
end