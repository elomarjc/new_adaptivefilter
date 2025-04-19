% Data for SNR values (before and after filtering for each experiment)
SNR_values = [
    -20.02, 3.25, 129.30, 113.32;   % Baseline on Sinus Wave
    -20.00, 5.82, 101.76, 91.62;    % Synthetic Heartbeat
    -0.27, 18.67, 16.69, 16.67;     % NHS with Hospital Noise
    -17.42, 23.43, 17.32, 17.29     % Your Own Experiment
];

% Define the filtering methods
filters = {'No Filtering', 'LMS', 'NLMS', 'RLS'};

% Define the experiment names
experiments = {'Baseline Sinus Wave', 'Synthetic Heartbeat', 'Company Data', 'Own Experiment'};

% Create a bar graph to compare SNR values
figure;
bar(SNR_values, 'grouped');
set(gca, 'xticklabel', experiments);
xlabel('Experiments');
ylabel('SNR (dB)');
%title('SNR Comparison Across Different Filtering Methods');
legend(filters, 'Location', 'best');
grid on;

% Adjust the graph for better visibility
ylim([-25 135]); % Adjust the Y-axis for better comparison
xtickangle(45); % Rotate x-axis labels for readability

% Crop the figure and save as PDF
tightfig();
saveas(gcf, 'SNR_comparison.pdf');

function tightfig()
    % Tighten the figure by removing excess whitespace
    set(gcf, 'Units', 'Inches');
    pos = get(gcf, 'Position');
    set(gcf, 'PaperUnits', 'Inches');
    set(gcf, 'PaperSize', [pos(3) pos(4)]);
    set(gcf, 'PaperPosition', [0 0 pos(3) pos(4)]);
end
