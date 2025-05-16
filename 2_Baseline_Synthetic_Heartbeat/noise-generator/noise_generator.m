% Parameters
fs = 8000;                 % Sample rate (Hz)
duration = 15;             % Duration (seconds)
f_mod = 2;                 % Square wave modulation frequency (Hz)
noise_bandwidth = 50;     % Bandwidth of the noise (Hz)
f_center = 500;           % Center frequency of the band
output_file = 'square_modulated_noise.wav';

% Time vector
t = 0:1/fs:duration-1/fs;

% Step 1: Generate square wave modulation
mod_signal = square(2*pi*f_mod*t);       % Values in [-1, 1]
mod_signal = (mod_signal + 1)/2;         % Now in [0, 1]

% Step 2: Generate white noise
white_noise = randn(size(t));            % Gaussian white noise

% Step 3: Apply bandpass filter to create narrowband noise
bpFilt = designfilt('bandpassiir', ...
    'FilterOrder', 8, ...
    'HalfPowerFrequency1', f_center - noise_bandwidth/2, ...
    'HalfPowerFrequency2', f_center + noise_bandwidth/2, ...
    'SampleRate', fs);
narrowband_noise = filter(bpFilt, white_noise);

% Step 4: Modulate noise with square wave
modulated_noise = mod_signal .* narrowband_noise;

% Normalize
modulated_noise = modulated_noise / max(abs(modulated_noise));

% Step 5: Save to WAV
audiowrite(output_file, modulated_noise, fs);
disp(['Square wave modulated noise saved to ', output_file]);
