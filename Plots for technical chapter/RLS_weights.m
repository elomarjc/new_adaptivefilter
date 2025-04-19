clear all

% Parameters
n = 50;  % Total number of points (arbitrary end point for demonstration)
i = 1:n;  % Time index
lambda_values = [0.5, 0.7, 0.9];  % Different lambda values

% Initialize figure
figure;
hold on;
title('Weighting Function \beta[n, i] for Different \lambda Values');
xlabel('i');
ylabel('\beta[n, i]');
grid on;

% Plot \beta[n, i] for each lambda value
for j = 1:length(lambda_values)
    lambda = lambda_values(j);
    beta = lambda .^ (n - i);  % Calculate weighting function
    plot(i, beta, 'LineWidth', 1.5, 'DisplayName', ['\lambda = ', num2str(lambda)]);
end

% Add legend
legend('show', 'Location', 'best');
hold off;


%%
clear all
close all
clc

% Parameters
n = 50;  % Total number of points (arbitrary end point for demonstration)
i = 1:n;  % Time index
lambda_values = [0.1, 0.7, 0.9];  % Different lambda values
colors = [0.85, 0.33, 0.1; 0, 0.5, 0; 0, 0.45, 0.74];  % Custom colors for each lambda value

% Initialize figure with specified size
figure;
set(gcf, 'Position', [1039 615 564 198]); % Set figure size in pixels (rounded values)

hold on;
xlabel('i', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('\beta[n, i]', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 10);  % Set axis font size
grid on;

% Plot \beta[n, i] for each lambda value
for j = 1:length(lambda_values)
    lambda = lambda_values(j);
    beta = lambda .^ (n - i);  % Calculate weighting function
    plot(i, beta, 'Color', colors(j, :), 'LineWidth', 2, 'DisplayName', ['\lambda = ', num2str(lambda)]);
end

% Customize axes and legend
ylim([-0.1 1.2]);  % Zoom out vertically to add margin
xlim([0 n + 5]);   % Zoom out horizontally to add margin

% Explicitly set legend with the correct entries only
legend({'\lambda = 0.1', '\lambda = 0.7', '\lambda = 0.9'}, 'Location', 'west', 'FontSize', 10);


% Additional styling for clarity
set(gca, 'XTick', [1, n], 'YTick', [0, 1]);  % Set discrete labels for x and y axes
set(gca, 'XTickLabel', {'1', 'n'}, 'YTickLabel', {'0', '1'});

% Plot dashed lines without adding them to the legend
plot([1, 1], [-0.1, 1], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');  % Dashed vertical line at i = 1
plot([n, n], [-0.1, 1], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');  % Dashed vertical line at i = n
plot([1, n], [1, 1], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');  % Dashed horizontal line at beta = 1

hold off;

% Crop the figure and save as PDF
tightfig();
saveas(gcf, 'RLS_weights.pdf');

function tightfig()
    % Tighten the figure by removing excess whitespace
    set(gcf, 'Units', 'Inches');
    pos = get(gcf, 'Position');
    set(gcf, 'PaperUnits', 'Inches');
    set(gcf, 'PaperSize', [pos(3) pos(4)]);
    set(gcf, 'PaperPosition', [0 0 pos(3) pos(4)]);
end