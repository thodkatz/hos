clc
clear
close all

% get sound
[sound,fs] = audioread('woman_o.wav');
sound = sound(:,1); % i don't know why I have two columns
[row,col] = size(sound); n = row;

estimatePitch = 0.004;

% start from 2sec for the frame
% assuming voice lasts over 2 seconds
Ts = 1/fs;
t = 0:Ts:(n-1)*Ts;
periods = 4;
startTime = 2;
startIndex = round(startTime/(Ts));
endIndex = round((startTime + periods*estimatePitch)/(Ts));
soundFrame = sound(startIndex:endIndex);
range = startIndex:endIndex;
nFrame = numel(soundFrame);
tshift = (-nFrame/2:nFrame/2-1)*(Ts)*1000;

figure
plot(t(1:numel(range))*1000,soundFrame);
xlabel('Time (ms)'),axis tight
title('Time signal')

% hamming Frame
soundHamming = soundFrame .* hamming(nFrame);
figure,plot(t(1:numel(range))*1000,soundHamming)
xlabel('Time (ms)'),title('Hamilton Frame applied'),axis tight

% cepstrum with hamming
soundCcepsHamming = fftshift(cceps(soundHamming));
plot(tshift,soundCcepsHamming),axis tight,title('Complex cepstrum with hamming')

% cpestrum without hamming
soundCcepsFrame = fftshift(cceps(soundFrame));
figure,plot(tshift,soundCcepsFrame),axis tight,title('Complex cepstrum without hamming')

% spot peak to estimate pitch
% real spectrum is a better tool for pitch estimation
soundRcepsHamming = fftshift(rceps(soundHamming));
figure,plot(tshift,soundRcepsHamming),axis tight,title('Real cepstrum with hamming')

% real spectrum without hamming
soundRcepsFrame = fftshift(rceps(soundFrame));
figure,plot(tshift,soundRcepsFrame),axis tight,title('Real cepstrum without hamming')

% spot pitch based on autocorr
acf = autocorr(soundHamming,numel(soundHamming)-1);
figure,plot(t(1:numel(range)),acf)
[pks,locs] = findpeaks(acf);
[max,index] = max(pks);
estimatePitch = locs(index)*Ts;

% deconvolution with hamming 
impulseResponse = icceps(ifftshift(soundCcepsHamming)) ./ hamming(nFrame);
figure,plot(t(1:numel(range)),impulseResponse),axis tight,title('Sanity check icceps time signal with hamming')

% deconvolution without hamming
impulseResponse = icceps(ifftshift(soundCcepsFrame));
figure,plot(t(1:numel(range)),impulseResponse),axis tight,title('Sanity check icceps time signal no hamming')

% design filter
filter = zeros(nFrame,1);
cutTime = estimatePitch;
cutIndex = round(cutTime/Ts);
centreIndex = round(numel(filter)/2);
centreRange = round(centreIndex-cutIndex:centreIndex+cutIndex);
filter(centreRange) = 1;

% liftering with hamming
soundCcepsHamming = filter .* soundCcepsHamming;
figure,plot(tshift,soundCcepsHamming),title('Liftering with hamming')
impulseResponse = (icceps(ifftshift(soundCcepsHamming))) ./ hamming(nFrame);
cutIndex = round(estimatePitch/Ts);
figure,plot(t(1:cutIndex)*1000,impulseResponse(1:cutIndex)),title('Impulse response after liftering hamming')

% liftering without hamming
% soundCcepsFrame = filter .* soundCcepsFrame;
% figure,plot(tshift,soundCcepsFrame),title('Liftering no hamming')
% impulseResponse = (icceps(ifftshift(soundCcepsFrame)));
% cutIndex = round(estimatePitch/Ts);
% figure,plot(t(1:cutIndex)*1000,impulseResponse(1:cutIndex)),title('Impulse response after liftering no hamming')

% reconstructing output
impulseResponse = impulseResponse(1:cutIndex);
impulseTrain = zeros(periods*numel(impulseResponse)+periods,1);
impulseTrain(1:numel(impulseResponse):end) = 1;
figure,plot(t(1:numel(impulseTrain)),impulseTrain)
soundFrameEst = conv(impulseResponse,impulseTrain);
figure,plot(t(1:numel(soundFrameEst))*1000,soundFrameEst);