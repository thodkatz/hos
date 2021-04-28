clear;
clc;
close all;
set(0,'defaultfigurecolor',[1 1 1])

N = 8192;
maxLag = 64;
maxLagPSD = 128;
NFFT = 2*maxLag+1;
M = 256;
samples = genSamples(N, 0);

% PSD
powerspectrum(samples, 0, maxLagPSD, 0);

% bispectrum L=64, M=256
hosAnalysis(samples, 256, maxLag, 0);

% repeat for different M = 512 and M = 128
hosAnalysis(samples, 128, maxLag, 0);
hosAnalysis(samples, 256, maxLag, 0);
hosAnalysis(samples, 512, maxLag, 0);

% repeat for numIter realizations and average the bispectrum and
% powerspectrum
numIter = 50;
psd = zeros(numIter, 2*maxLagPSD+1);
bispPa = zeros(numIter, NFFT, NFFT);
bispNo = zeros(numIter, NFFT, NFFT);
bispDirect = zeros(numIter, M, M);
for k = 1:numIter
    samples = genSamples(N, 0);
    [psd(k,:),bispPa(k,:,:), bispNo(k,:,:), frequency, bispDirect(k,:,:), ...
        waxis] = hosAnalysis(samples, M, maxLag, 0);
end

% calculate means
psdMean = mean(psd, 1);
bispPaMean = mean(abs(bispPa), 1);
bispNoMean = mean(abs(bispNo), 1);
bispDirectMean = mean(abs(bispDirect), 1);
bispPaMean = reshape(bispPaMean, NFFT, NFFT);
bispNoMean = reshape(bispNoMean, NFFT, NFFT);
bispDirectMean = reshape(bispDirectMean, M, M);

% calculate variance coefficient
bispPaStd = std(abs(bispPa), 0, 1);
bispNoStd = std(abs(bispNo), 0, 1);
bispDirectStd = std(abs(bispDirect), 0, 1);
bispPaStd = reshape(bispPaStd, NFFT, NFFT);
bispNoStd = reshape(bispNoStd, NFFT, NFFT);
bispDirectStd = reshape(bispDirectStd, M, M);

% PLOTS
display = 1;
if display ~= 0
    fs = 1;
    fshift = (-maxLagPSD:maxLagPSD)*(fs/(2*maxLagPSD+1)); % zero-centered frequency range
    
    figure();
    plot(fshift, psdMean);
    title("Average PSD");
    
    figure();
    plot(std(psd, 0, 1));
    title("Average PSD standard deviation");
end

if display ~= 0
    figure();
    subplot(211);
    hold on;
    plot(frequency(maxLag+1:end), frequency(maxLag+1:end), 'color', 'red');
    contour(frequency(maxLag+1:end), frequency(maxLag+1:end), ...
        bispPaMean(maxLag+1:end, maxLag+1:end), 8);
    title("Average Bispectrum Indirect Parzen");
    subplot(212);
    mesh(frequency(maxLag+1:end), frequency(maxLag+1:end), ...
        bispPaMean(maxLag+1:end,maxLag+1:end)), colorbar;
    
    figure(); 
    hold on;
    plot(frequency(maxLag+1:end), frequency(maxLag+1:end), 'color', 'red');
    contour(frequency(maxLag+1:end), frequency(maxLag+1:end), ...
        bispPaStd(maxLag+1:end,maxLag+1:end),8), colorbar;
    title("Average Bispectrum Indirect Parzen standard deviation");
    
    figure();
    subplot(211);
    hold on;
    plot(frequency(maxLag+1:end), frequency(maxLag+1:end), 'color', 'red');
    contour(frequency(maxLag+1:end), frequency(maxLag+1:end), ...
        bispNoMean(maxLag+1:end, maxLag+1:end), 8), colorbar;
    title("Average Bispectrum Indirect No window");
    subplot(212);
    mesh(frequency(maxLag+1:end), frequency(maxLag+1:end), ...
        bispNoMean(maxLag+1:end, maxLag+1:end)), colorbar;
    
    figure(); 
    hold on;
    plot(frequency(maxLag+1:end), frequency(maxLag+1:end), 'color', 'red');
    contour(frequency(maxLag+1:end), frequency(maxLag+1:end), ...
        bispNoStd(maxLag+1:end,maxLag+1:end),8), colorbar;
    title("Average Bispectrum Indirect No window standard deviation");
    
    figure();
    subplot(211);
    hold on;
    plot(waxis(M/2:end), waxis(M/2:end), 'color', 'red');
    contour(waxis(M/2:end), waxis(M/2:end), ...
        abs(bispDirectMean(M/2:end, M/2:end)), 8);
    title("Average Bispectrum Direct");
    subplot(212);
    mesh(waxis(M/2:end), waxis(M/2:end), ...
        bispDirectMean(M/2:end, M/2:end)), colorbar;
    
    figure(); 
    hold on;
    plot(waxis(M/2:end), waxis(M/2:end), 'color', 'red');
    contour(waxis(M/2:end), waxis(M/2:end), ...
        bispDirectStd(M/2:end,M/2:end),8), colorbar;
    title("Average Bispectrum Direct standard deviation");
end

function samples = genSamples(N, display)

% GENSAMPLES Generate samples for the given time series
% display ~= 0 Display time series sample
% N The number of samples

    % frequency armonics
    lambda = zeros(1, 6);
    lambda(1) = 0.12;
    lambda(2) = 0.3;
    lambda(3) = lambda(1) + lambda(2);
    lambda(4) = 0.19;
    lambda(5) = 0.17;
    lambda(6) = lambda(4) + lambda(5);

    % phase armonics
    phi = zeros(1, 6);
    phi(1) = 2*pi*rand(1, 1);
    phi(2) = 2*pi*rand(1, 1);
    phi(3) = phi(1) + phi(2);
    phi(4) = 2*pi*rand(1, 1);
    phi(5) = 2*pi*rand(1, 1);
    phi(6) = phi(4) + phi(5);

    % create random signal
    syms t;
    x = cos(2*pi*lambda(1)*t + phi(1));
    for i = 2:6
        x = x + cos(2*pi*lambda(i)*t + phi(i));
    end
    x = matlabFunction(x);
    
    if display ~= 0
        figure();
        plot(100:1:200, x(100:1:200));
        title("Time series");
        xlabel("k");
    end
    
    samples = x(0:1:N-1);
end

function primary = primaryArea(array2d)

% Get a 2d array and keep only the lower half

[row, col] = size(array2d);

filter = ones(row, col);
filter = tril(filter)';
%filter = flip(filter);

primary = array2d .* filter;

end