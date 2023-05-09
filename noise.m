% Sampling frequency
Fs = 2000;

% Signal duration (in seconds)
duration = 5;

% Frequency range (in Hz)
fmin = 180;
fmax = 220;

% Number of frequency components
N = 50;

% Generate a vector of N equally spaced frequencies within the specified range
frequencies = linspace(fmin, fmax, N);

% Generate random amplitudes for each frequency component
amplitudes = randn(N, 1);

% Generate a time vector
t = linspace(0, duration, Fs * duration);

% Initialize the signal
signal = zeros(1, length(t));

% Loop over each frequency component and add it to the signal
for i = 1:N
    % Generate a sinusoidal signal with the current frequency and random amplitude
    component = amplitudes(i) * sin(2 * pi * frequencies(i) * t);
    
    % Add the component to the signal
    signal = signal + component;
end

% Normalize the signal to have unit power
signal = signal / sqrt(sum(signal.^2));

% Plot the signal
figure(1)
plot(t, signal);
xlabel('Time (s)');
ylabel('Amplitude');
title('Band-limited white noise');