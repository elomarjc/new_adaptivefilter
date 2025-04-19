%%
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
mus = 1e-3 * (2.^(0:9));  % Step sizes for NLMS
ps = 1:10;  % Filter orders

for folder = folders'
    if folder.isdir
        folderPath = fullfile(dataPath, folder.name);
        noisyFile = fullfile(folderPath, "primary.wav");
        referenceFile = fullfile(folderPath, "secondary.wav");
        cleanFile = fullfile(folderPath, "ZCH0019.wav");
        
        if isfile(noisyFile) && isfile(referenceFile) && isfile(cleanFile)
            [d, Fs] = audioread(noisyFile);
            [x, ~] = audioread(referenceFile);
            [clean, ~] = audioread(cleanFile);

            % Loop through different configurations for NLMS
            SNRs = zeros(length(mus), length(ps));
            for i = 1:length(mus)
                for j = 1:length(ps)
                    
                    % Parameters
                    mu = mus(i);
                    p = ps(j);
                    
                    % Initialization for NLMS
                    w = zeros(p, 1);
                    x_n = zeros(p, 1);  % Signal fragment to be filtered
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
                    
                    % Track best SNR and parameters
                    if snr_NLMS > bestSNR_NLMS
                        bestSNR_NLMS = snr_NLMS;
                        bestParams_NLMS = [p, mu];
                    end
                end
            end

            % Display the best parameters
            fprintf("Best NLMS: Filter Order = %d, Step Size = %.4f, SNR = %.2f dB\n", bestParams_NLMS(1), bestParams_NLMS(2), bestSNR_NLMS);

            % Plot time evolution of the filter weights for the best SNR configuration
            opt_mu = bestParams_NLMS(2);  % Optimal step-size parameter
            opt_p = bestParams_NLMS(1);  % Optimal filter order
            
            % Initialization for optimal configuration
            w = zeros(opt_p, length(d) + 1);
            x_n = zeros(opt_p, 1);  % Signal fragment to be filtered
            e = zeros(length(d), 1);

            % Filtering with optimal configuration
            for k = 1:length(d)
                x_n(2:end) = x_n(1:end-1);
                x_n(1) = x(k);
                sig = sum(x_n.^2);
                e(k) = d(k) - w(:, k)' * x_n;
                w(:, k + 1) = w(:, k) + (opt_mu / (sig + 1e-10)) * x_n * e(k);
            end

            % Plot weights evolution
            figure;
            plot(w');
            grid on;
            title('NLMS Weights Evolution');
            xlabel('Cycle (n)');
            ylabel('Magnitude');
            legend(arrayfun(@(i) sprintf('w_%d(n)', i - 1), 1:opt_p, 'UniformOutput', false));
            
            % Plot signal and filtered output
            figure;
            plot(clean);
            hold on;
            plot(e, 'r--');
            hold off;
            legend('Clean', 'Filtered');
            grid on;
            xlabel('Time (n)');
            ylabel('Amplitude');
        end
    end
end



function tightfig()
    % Tighten the figure by removing excess whitespace
    set(gcf, 'Units', 'Inches');
    pos = get(gcf, 'Position');
    set(gcf, 'PaperUnits', 'Inches');
    set(gcf, 'PaperSize', [pos(3) pos(4)]);
    set(gcf, 'PaperPosition', [0 0 pos(3) pos(4)]);
end