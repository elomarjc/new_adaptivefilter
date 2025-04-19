%% Clear workspace and close figures
clear all;
close all;
clc;

%% Load the primary and secondary signals (noise and reference)
[primary, fs] = audioread("C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Matlab\adaptivefilterMT\3_Baseline_Company_Data\Experiment_Data_vs2\Hospital Ambient Noises\NHS\1\primary.wav");
[noise, ~] = audioread("C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Matlab\adaptivefilterMT\3_Baseline_Company_Data\Experiment_Data_vs2\Hospital Ambient Noises\NHS\1\secondary.wav");
%[clean, ~] = audioread("C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Data from AI heathway\Data_ANC\Experiment_Data\Hospital Ambient Noises\NHS\1\ZCH0019.wav");
% Clean signal is not needed for ICA, only for SNR calculation after filtering

% Ensure the signals are of the same length
minLength = min([length(primary), length(noise)]);  % Only use primary and noise for ICA
primary = primary(1:minLength);
noise = noise(1:minLength);

%% Center the data (as required by PowerICA)
Y = [primary'; noise'];  % Create a 2-by-n mixture matrix
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

%% Calculate SNR (optional)
% Assuming you have a clean signal available for SNR calculation:
% filtered_signal_1 = S_est(1, :); % Separated first source (heart/lung sound)
% filtered_signal_2 = S_est(2, :); % Separated second source (noise)

% SNR calculation after ICA separation (use the clean signal for SNR)
% filtered_SNR_1 = 10 * log10(sum(clean.^2) / sum((clean - filtered_signal_1).^2));
% filtered_SNR_2 = 10 * log10(sum(clean.^2) / sum((clean - filtered_signal_2).^2));

%% Export separated signals after ICA (optional)
audiowrite('Separated_Signal_1.wav', S_est(1, :), fs);  % Separated source 1 (heart/lung)
audiowrite('Separated_Signal_2.wav', S_est(2, :), fs);  % Separated source 2 (noise)

fprintf('Separated signals saved as audio files.\n\n');

%% Visualize Results
figure;
subplot(2, 1, 1);
plot(primary);
title('Original Noisy Signal');
xlabel('Sample Number');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(S_est(1, :));
title('Estimated Source Signal 1');
xlabel('Sample Number');
ylabel('Amplitude');

% Save the figure as a PDF
saveas(gcf, 'PowerICA_Separation_Results.pdf');
