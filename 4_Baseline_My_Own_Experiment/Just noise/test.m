clear all
close all

[u, fs] = audioread("C:\Users\eloma\Desktop\new_adaptivefilter\4_Baseline_My_Own_Experiment\Just noise\Audio 1.wav");   %noise + clean signal
[d, ~] = audioread("C:\Users\eloma\Desktop\new_adaptivefilter\4_Baseline_My_Own_Experiment\Just noise\Audio 2.wav");


error_signal = u - d;
plot(error_signal)
xlabel('Sample Number')
ylabel('Amplitude')
title('Error Signal (Mic 1 - Mic 2)')
grid on

% Set figure size and position
set(gcf, 'Units', 'inches', 'Position', [7.739583333333333,8.675,5.833333333333333,1.210416666666667]);

% Save the figure
tightfig();

% Construct filename
saveas(gcf, 'Error signal of 2 mics.pdf');

%% Save Figure trimmed
function tightfig()
    set(gcf, 'Units', 'Inches');
    pos = get(gcf, 'Position');
    set(gcf, 'PaperUnits', 'Inches');
    set(gcf, 'PaperSize', [pos(3) pos(4)]);
    set(gcf, 'PaperPosition', [0 0 pos(3) pos(4)]);
end