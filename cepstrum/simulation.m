clc
clear
close all

% get sound
[sound,fs] = audioread('woman_o.wav');
sound = sound(:,1); % i don't know why I have two columns
[row,col] = size(sound); n = row;

estimatePitch = 0.005;

% start from 2sec for the frame
% assuming voice lasts over 2 seconds
Ts = 1/fs;
t = 0:Ts:(n-1)*Ts;
periods = 5;
startTime = 2;
startIndex = round(startTime/(Ts));
endIndex = round((startTime + periods*estimatePitch)/(Ts));
soundFrame = sound(startIndex:endIndex);
range = startIndex:endIndex;
nFrame = numel(soundFrame);

figure
subplot(2, 2, 1)
plot(t(1:numel(range))*1000,soundFrame);
xlabel('Time (ms)'),axis tight
title('Time signal')

% hamming Frame
subplot(2, 2, 3)
soundFrame = soundFrame .* hamming(nFrame);
plot(t(1:numel(range))*1000,soundFrame)
xlabel('Time (ms)'),title('Hamilton Frame applied'),axis tight

% ceps and reps in full signal
subplot(2, 2, 2)
rcep = fftshift(rceps(sound));
tshift = (-nFrame/2:nFrame/2-1)*(Ts)*1000;
centreIndex = round(numel(rcep)/2);
centreRange = round(centreIndex-nFrame/2:centreIndex+nFrame/2-1);
rcep = rcep(centreRange);
plot(tshift,rcep)
xlabel('quefrency (ms)'),title('Real cepstrum')

subplot(2, 2, 4)
ccep = fftshift(cceps(sound));
tshift = (-nFrame/2:nFrame/2-1)*(Ts)*1000;
centreIndex = round(numel(ccep)/2);
centreRange = round(centreIndex-nFrame/2:centreIndex+nFrame/2-1);
ccep = ccep(centreRange);
plot(tshift,ccep)
xlabel('quefrency (ms)'),title('Complex cepstrum')

% FFT
Y = fftshift(fft(soundFrame));
fshift = (-nFrame/2:nFrame/2-1)*(fs/nFrame);
figure,plot(fshift,abs(Y)/nFrame)
xlabel('frequency (Hz)'),title('FFT')

% design filter
filter = zeros(nFrame,1);
cutTime = estimatePitch;
cutIndex = round(cutTime/Ts);
centreIndex = round(numel(filter)/2);
centreRange = round(centreIndex-cutIndex:centreIndex+cutIndex);
filter(centreRange) = 1;

% rceps in frame 
rcepFrame = fftshift(rceps(soundFrame));
tshift = (-nFrame/2:nFrame/2-1)*(Ts)*1000;
figure,plot(tshift,rcepFrame)
xlabel('quefrency (ms)'),title('Real cepstrum Frame')

% cceps in frame
ccepFrame = fftshift(cceps(soundFrame));
tshift = (-nFrame/2:nFrame/2-1)*(Ts)*1000;
figure
hold on
plot(tshift,filter);
plot(tshift,ccepFrame)
xlabel('quefrency (ms)'),title('Complex cepstrum Frame')

% liftering ccep in frame
impulseCcepFrame = filter .* ccepFrame;
inputCcepFrame = ccepFrame - impulseCcepFrame;
figure
hold on
yline(6) % just for scaling
yline(-6)
plot(tshift,inputCcepFrame),title('Liftering impulse train')
figure
plot(tshift,impulseCcepFrame),title('Liftering impulse response')

% estimate impulse ccep frame signal
impulseHat = icceps(ifftshift(impulseCcepFrame));
figure,plot(impulseHat),title('Impulse response estimation')
inputHat = icceps(ifftshift(inputCcepFrame));
figure,plot(inputHat),title('Impulse train estimation')

sound = conv(inputHat,impulseHat,'same');
figure,plot(t(1:numel(range))*1000,sound),title('Time signal estimation')
