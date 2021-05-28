clc
clear
close all

% get sound
[sound,fs] = audioread('a_good_me.wav');
sound = sound(:,1); % i don't know why I have two columns
[row,col] = size(sound); N = row;

% f0 = pitch(sound,fs,'Method','CEP');
% figure
% plot(f0)

estimatePitch = 0.005;

% start from 2sec for the frame
% assuming voice lasts over 2 seconds
Ts = 1/fs;
t = 0:Ts:(N-1)*Ts;
periods = 4;
tStart   = round(2/(Ts));
tEnd     = round((2 + periods*estimatePitch)/(Ts));
soundWindow = sound(tStart:tEnd);
range     = tStart:tEnd;
N = numel(soundWindow);

figure
subplot(2, 2, 1)
plot(t(1:numel(range))*1000,soundWindow);
xlabel('Time (ms)'),axis tight
title('Time signal')

% hamming window
subplot(2, 2, 3)
soundWindow = soundWindow .* hamming(N);
plot(t(1:numel(range))*1000,soundWindow)
xlabel('Time (ms)'),title('Hamilton window applied'),axis tight

% ceps and reps in full signal
halfCeps = round(numel(range)/2);
subplot(2, 2, 2)
rcep = rceps(sound);
plot(t(1:halfCeps)*1000,rcep(1:halfCeps))
xlabel('quefrency (ms)'),title('Real cepstrum')

subplot(2, 2, 4)
rcep = cceps(sound);
plot(t(1:halfCeps)*1000,rcep(1:halfCeps))
xlabel('quefrency (ms)'),title('Complex cepstrum')

figure
Y = fftshift(fft(soundWindow));
fshift = (-N/2:N/2-1)*(fs/N);
plot(fshift,abs(Y)/N)
xlabel('frequency (Hz)'),title('FFT')

% ceps and reps in frame 
figure
subplot(2, 1, 1)
rcep = rceps(soundWindow);
plot(t(1:halfCeps)*1000,rcep(1:halfCeps))
xlabel('quefrency (ms)'),title('Real cepstrum')

subplot(2, 1, 2)
rcep = cceps(soundWindow);
plot(t(1:halfCeps)*1000,rcep(1:halfCeps))
xlabel('quefrency (ms)'),title('Complex cepstrum')