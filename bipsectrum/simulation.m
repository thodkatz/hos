clear;
clc;
close all;

% frequency armonics
lambda = zeros(1, 6);
lambda(1) = 0.12;
lambda(2) = 0.3;
lambda(3) = lambda(1) + lambda(2);
lambda(4) = 0.19;
lambda(5) = 0.17;
lambda(6) = lambda(4) + lambda(5);
lambdas = lambda'

% phase armonics
phi = zeros(1, 6);
phi(1) = 2*pi*rand(1, 1);
phi(2) = 2*pi*rand(1, 1);
phi(3) = phi(1) + phi(2);
phi(4) = 2*pi*rand(1, 1);
phi(5) = 2*pi*rand(1, 1);
phi(6) = phi(4) + phi(5);

% samples
N = 2^13;

% ensure uniform distribution
samples = zeros(1, N);
for i = 1:N
    random = 2*pi*rand(1, 1);
    samples(i) = random;
end
%histogram(samples);

% create random signal
syms k;
x = cos(2*pi*lambda(1)*k + phi(1));
for i = 2:6
    x = x + cos(2*pi*lambda(i)*k + phi(i));
end
x = matlabFunction(x);
%figure();
%plot(1:1:N, x(0:1:N-1));
%title("Time series");

samples = x(0:1:N-1);

% autocorrelation hosa toolbox
L2 = 128;
cum2 = cumest(samples, 2, L2, 256, 0, 'unbiased');

% normalizing with the variance
% should the autocorrelation be normalized for power spectrum?
%cum2 = cum2 ./ var(samples);

figure();
%plot(-L2:L2, cum2);
plot(0:L2, cum2(L2+1:end)); % symmetric
title("Autocorrelation cumest");

% autocorrelation via xcorr
acf = xcorr(samples, L2, 'unbiased');
figure();
%plot(-L2:L2, acf);
plot(0:L2, acf(L2+1:end));
title("Autocorrelation xcorr");

% power spectrum
fs = 1;
fshift = (-L2:L2)*(fs/(2*L2+1)); % zero-centered frequency range
powershift = abs(fftshift(fft(acf)));
figure();
%plot(fshift,powershift);
plot(fshift(L2+1:end), powershift(L2+1:end));
title("Power spectral density autocorr");

powershift = abs(fftshift(fft(cum2)));
figure();
plot(fshift, powershift);
title("Power spectral density cumest 1");

% test for fft. Does it matter getting only the half positive vs all
%powershift = abs(fftshift(fft(cum2(L2:2*L2-1))));
%figure();
%plot(-L2/2:L2/2-1, powershift);
%title("Power spectral density cumest 2");

% third order cumulant
L2 = 128;
M = 256;
for k = -L2:L2
cum3(:,k+L2+1) = cumest(samples, 3, L2, 0, 0, 'unbiased', k); % won't work for M as sample length
end

% should the autocorrelation be normalized with kyrtosis?
%kyrtosis = cum3(L2+1,L2+1) % c3(0,0)
%cum3 = cum3 ./ kyrtosis;

figure();
%subplot(211);
mesh(-L2:L2, -L2:L2, cum3);
title("Third order cumulant cumest");
%subplot(212);
%contour(-L2:L2, -L2:L2, cum3, 8);

% bispectrum indirect
NFFT = 2*L2 + 1;
%NFFT = L2; % this creates a better contour tho, why?
bispeci(samples, L2, M, 0, 'unbiased', NFFT, -1);
bispeci(samples, L2, M, 0, 'unbiased', NFFT, 0);

% bispectrum alternative
figure();
[bisp, frequency, cum3alt] = bisp3cum(samples,1,L2,'pa', 'unbiased');
figure();
mesh(-L2:L2, -L2:L2, cum3alt);
title("Third order cumulant bisp3cum");
figure();
subplot(211);
contour(frequency, frequency, abs(bisp));
title("Bispectrum bisp3cum parzen");
subplot(212);
mesh(frequency(L2+1:end), frequency(L2+1:end), abs(bisp(L2+1:end,L2+1:end)));

figure();
[bisp, frequency, cum3alt] = bisp3cum(samples,1,L2,'none', 'unbiased');
figure();
subplot(211);
contour(frequency, frequency, abs(bisp));
title("Bispectrum bisp3cum no window");
subplot(212);
mesh(frequency(L2+1:end), frequency(L2+1:end), abs(bisp(L2+1:end, L2+1:end)));

% bispectrum direct
%NFFT = 
%wind = ones(); % rectangular window - no window
%bispecd(samples)