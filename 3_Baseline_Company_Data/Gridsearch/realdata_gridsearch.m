clc;
clear;
close all;

% Define paths
dataPath = "C:\Users\eloma\Desktop\Universitet\OneDrive - Aalborg Universitet\Universitet\9. Semester - ES9\Long Thesis\Data from AI heathway\Data_ANC\Experiment_Data\Hospital Ambient Noises\NHS\1";
folders = dir(fullfile(dataPath, "*"));

% Define SNR improvement tracking for NLMS
bestSNR_NLMS = -Inf;
bestParams_NLMS = [];

% Define NLMS parameter search space
mus = 1e-4 * (2.^(0:9));  % Step sizes for NLMS
ps = 1:25;  % Filter lengths

for folder = folders'
    if folder.isdir
        folderPath = fullfile(dataPath, folder.name);
        noisyFile = fullfile(folderPath, "primary.wav");
        referenceFile = fullfile(folderPath, "secondary.wav");
        cleanFile = fullfile(folderPath, "ZCH0019.wav");

        if isfile(noisyFile) && isfile(referenceFile) && isfile(cleanFile)
            [x, Fs] = audioread(noisyFile);
            [d, ~] = audioread(referenceFile);
            [clean, ~] = audioread(cleanFile);

            % Loop through different configurations for NLMS
            SNRs = zeros(length(mus), length(ps));
            for i = 1:length(mus)
                for j = 1:length(ps)
                    mu = mus(i);
                    p = ps(j);

                    % Initialization
                    w = zeros(p, 1);
                    x_n = zeros(p, 1);
                    e = zeros(length(d), 1);

                    % NLMS Filtering
                    for k = 1:length(d)
                        x_n(2:end) = x_n(1:end-1);
                        x_n(1) = x(k);
                        sig = sum(x_n.^2);
                        e(k) = d(k) - w' * x_n;
                        w = w + (mu / (sig + 1e-10)) * x_n * e(k);
                    end

                    % SNR computation
                    snr_NLMS = 10 * log10(sum(clean.^2) / sum((clean - e).^2));
                    SNRs(i, j) = snr_NLMS;

                    % Track best configuration
                    if snr_NLMS > bestSNR_NLMS
                        bestSNR_NLMS = snr_NLMS;
                        bestParams_NLMS = [p, mu];
                    end
                end
            end

            % Display best parameters
            fprintf("Best NLMS: Filter Order = %d, Step Size = %.4f, SNR = %.2f dB\n", bestParams_NLMS(1), bestParams_NLMS(2), bestSNR_NLMS);

            % Re-run NLMS with best configuration
            opt_mu = bestParams_NLMS(2);
            opt_p = bestParams_NLMS(1);
            w = zeros(opt_p, length(d) + 1);
            x_n = zeros(opt_p, 1);
            e = zeros(length(d), 1);

            for k = 1:length(d)
                x_n(2:end) = x_n(1:end-1);
                x_n(1) = x(k);
                sig = sum(x_n.^2);
                e(k) = d(k) - w(:, k)' * x_n;
                w(:, k + 1) = w(:, k) + (opt_mu / (sig + 1e-10)) * x_n * e(k);
            end

            % === PLOTS ===

            % 1. Weight Evolution
            figure;
            plot(w');
            grid on;
            title('NLMS Weights Evolution');
            xlabel('Cycle (n)');
            ylabel('Magnitude');
            legend(arrayfun(@(i) sprintf('w_%d(n)', i - 1), 1:opt_p, 'UniformOutput', false));

            % 2. Time-domain Signals: Clean vs Filtered
            figure;
            plot(clean);
            hold on;
            plot(e, 'r--');
            hold off;
            legend('Clean', 'Filtered');
            grid on;
            xlabel('Time (n)');
            ylabel('Amplitude');
            title('Clean vs Filtered Output');

            % 3. Best SNR per Filter Length
            bestSNR_per_p = zeros(1, length(ps));
            for j = 1:length(ps)
                bestSNR_per_p(j) = max(SNRs(:, j));
            end

            figure;
            plot(ps, bestSNR_per_p,'-o', 'LineWidth', 1.5, 'DisplayName', 'NLMS');
            grid on;
            xlabel('Filter Length (M)');
            ylabel('Best SNR (dB)');
            title('Best SNR vs Filter Length (NLMS)');

            % 4. Frequency Response of Final Weights
            n_fft = 257;  % Number of FFT points
            x_range = linspace(0, Fs/2, n_fft);  % Frequency axis from 0 to fs/2
            
            wf = w(:, end);
            fr = 20 * log10(abs(freqz(wf, 1, n_fft)));
            figure;
            plot(x_range, fr);
            xlabel('Frequency (Hz)');
            ylabel('Magnitude (dB)');
            title('Frequency Response of NLMS Filter');
            grid on;

            % 5. Mel Spectrograms
            figure;
            subplot(2,1,1);
            melSpectrogram(clean, Fs);
            title('Mel Spectrogram - Clean Signal');

            subplot(2,1,2);
            melSpectrogram(e, Fs);
            title('Mel Spectrogram - NLMS Filtered Signal');

            % === Save Output Audio ===
            e = e / max(abs(e));
            outputFile = "Filtered_NLMS.wav";
            audiowrite(outputFile, e, Fs);
            fprintf("Filtered audio saved: %s\n", outputFile);
        end
    end
end

function tightfig()
    set(gcf, 'Units', 'Inches');
    pos = get(gcf, 'Position');
    set(gcf, 'PaperUnits', 'Inches');
    set(gcf, 'PaperSize', [pos(3) pos(4)]);
    set(gcf, 'PaperPosition', [0 0 pos(3) pos(4)]);
end
