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

% get time interval 2sec - 2.002 sec
% assuming voice lasts over 2 seconds
Ts = 1/fs;
t = 0:Ts:(N-1)*Ts;
periods = 20;
tStart   = round(2/(Ts));
tEnd     = round((2 + periods*estimatePitch)/(Ts));
tEndHalf = round((2 + periods*estimatePitch/2)/(Ts));
soundWindow = sound(tStart:tEnd);
N = numel(soundWindow);

figure
subplot(2, 2, 1)
plot(t(tStart:tEnd)*1000,soundWindow);
xlabel('Time (ms)'),axis tight
title('Time signal')

% hamming window
subplot(2, 2, 3)
soundWindow = soundWindow .* hamming(N);
plot(t(tStart:tEnd),soundWindow)
title('Windowed signal')
axis tight

subplot(2, 2, 2)
c = cceps(soundWindow);
range = tStart:tEndHalf;
plot(t(range)*1000,c(1:numel(range)))
xlabel('Time (ms)'),axis tight
title('Complex cepstrum')

subplot(2, 2, 4)
r = rceps(soundWindow);
plot(t(range)*1000,r(1:numel(range)))
xlabel('Time (ms)'),axis tight
title('Real cepstrum')

figure
Y = fftshift(fft(soundWindow));
fshift = (-N/2:N/2-1)*(fs/N);
plot(fshift,abs(Y)/N)

figure
rcep = rceps(sound);
plot(t(1:numel(range)),rcep(1:numel(range)))
title('Real cepstrum all signal Hamming window')

figure
plot(rceps(sound))

figure
plot(cceps(sound))