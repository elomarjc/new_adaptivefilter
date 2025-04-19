close all

% Parameters
mu_values = [0.1, 0.2, 0.5, 2]; % Step sizes including mu = 2
lambda = 1;
n = 0:50; % Number of adaptation cycles

% Transparency level for mu = 2
alpha = 0.3;

% Create figure
figure('Position', [100, 100, 600, 250]);
hold on;
grid on;

% Plot for each mu
for i = 1:length(mu_values)
    mu = mu_values(i);
    y = (1 - mu * lambda) .^ n; % Compute values
    
    if mu == 2
        % Plot for mu = 2 with grey color and transparency
        h = plot(n, y, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5, 'DisplayName', sprintf('\\mu = %.1f', mu));
        h.Color(4) = alpha; % Apply transparency to the greyed-out plot
    else
        % Plot for the other mu values with default color (no color specified)
        plot(n, y, 'LineWidth', 1.5, 'DisplayName', sprintf('\\mu = %.1f', mu));
    end
end

% Customize plot
xlabel('Number of adaptation cycles, n');
ylabel('$(1 - \mu_k \lambda)^n$', 'FontSize', 14, 'Interpreter', 'latex');
%title('Convergence behavior for different step sizes \mu');
legend show;
ylim([-1.5 1.5]); % Adjusted for oscillation visibility

% Export figure to PDF without cropping
set(gca, 'LooseInset', [0.1, 0, 0, 0]); % Small padding to prevent clipping
exportgraphics(gcf, 'LMS_mu_behavior.pdf', 'ContentType', 'vector', ...
    'BackgroundColor', 'none', 'Resolution', 1000);
