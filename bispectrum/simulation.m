clear;
clc;
close all;

display = 0;
N = 8192;
samples = genSamples(N, display);

maxLag = 128;
M = 256;

rng(0);

% repeat for different M = 512 and M = 128
hosAnalysis(samples, 128, maxLag, 0);
hosAnalysis(samples, 256, maxLag, 0);
hosAnalysis(samples, 512, maxLag, 0);

% repeat for 50 realizations and get the mean
numIter = 4;
psd = zeros(numIter, 2*maxLag + 1);
bispPa = zeros(numIter, 2*maxLag + 1, 2*maxLag + 1);
bispNo = zeros(numIter, 2*maxLag + 1, 2*maxLag + 1);
bispDirect = zeros(numIter, M, M);
for k = 1:numIter
    samples = genSamples(N, display);
    [psd(k,:),bispPa(k,:,:), bispNo(k,:,:), frequency, bispDirect(k,:,:), waxis] = ...
        hosAnalysis(samples, M, maxLag, 1);
end

bispPaStd = std(abs(bispPa), 0, 1)./mean(abs(bispPa));
bispNoStd = std(abs(bispNo), 0, 1)./ mean(abs(bispNo));
bispDirectStd = std(abs(bispDirect), 0, 1) ./ mean(abs(bispDirect));
bispPaStd = reshape(bispPaStd, 2*maxLag+1, 2*maxLag+1);
bispNoStd = reshape(bispNoStd, 2*maxLag+1, 2*maxLag+1);
bispDirectStd = reshape(bispDirectStd, M, M);

psdMean = mean(psd, 1);
bispPaMean = mean(abs(bispPa), 1);
bispNoMean = mean(abs(bispNo), 1);
bispDirectMean = mean(abs(bispDirect), 1);
bispPaMean = reshape(bispPaMean, 2*maxLag+1, 2*maxLag+1);
bispNoMean = reshape(bispNoMean, 2*maxLag+1, 2*maxLag+1);
bispDirectMean = reshape(bispDirectMean, M, M);

%%% PLOTS
display = 1;
if display ~= 0
    fs = 1;
    fshift = (-maxLag:maxLag)*(fs/(2*maxLag+1)); % zero-centered frequency range
    
    figure();
    plot(fshift, psdMean);
    title("Power spectral density cumest");
    
    figure();
    plot(var(psd, 0, 1));
    title("PSD standard deviation");
end

if display ~= 0
    figure();
    subplot(211);
    contour(frequency, frequency, bispPaMean, 8);
    title("Bispectrum bisp3cum parzen");
    subplot(212);
    mesh(frequency(maxLag+1:end), frequency(maxLag+1:end), ...
        bispPaMean(maxLag+1:end,maxLag+1:end));
        
    figure(); 
    contour(frequency(maxLag+1:end), frequency(maxLag+1:end), ...
        bispPaStd(maxLag+1:end,maxLag+1:end),8);
    title("Bispectrum Parzen coefficient of variation");
    
    figure();
    subplot(211);
    contour(frequency, frequency, bispNoMean, 8);
    title("Bispectrum bisp3cum no window");
    subplot(212);
    mesh(frequency(maxLag+1:end), frequency(maxLag+1:end), ...
        bispNoMean(maxLag+1:end, maxLag+1:end));
    
    figure(); 
    contour(frequency(maxLag+1:end), frequency(maxLag+1:end), ...
        bispNoStd(maxLag+1:end,maxLag+1:end),8);
    title("Bispectrum No window coefficient of variation");
   
    figure();
    subplot(211);
    contour(waxis, waxis, abs(bispDirectMean), 8);
    title("Bispectrum bisp3cum Direct");
    subplot(212);
    mesh(waxis(M/2:end), waxis(M/2:end), ...
        bispDirectMean(M/2:end, M/2:end));
    
    figure(); 
    contour(waxis(M/2:end), waxis(M/2:end), ...
        bispDirectStd(M/2:end,M/2:end),8);
    title("Bispectrum Direct coefficient of variation");
end

% if display ~= 0 then display time series sample
function samples = genSamples(N, display)
    % frequency armonics
    lambda = zeros(1, 6);
    lambda(1) = 0.12;
    lambda(2) = 0.3;
    lambda(3) = lambda(1) + lambda(2);
    lambda(4) = 0.19;
    lambda(5) = 0.17;
    lambda(6) = lambda(4) + lambda(5);
    lambdas = lambda';

    % phase armonics
    phi = zeros(1, 6);
    phi(1) = 2*pi*rand(1, 1);
    phi(2) = 2*pi*rand(1, 1);
    phi(3) = phi(1) + phi(2);
    phi(4) = 2*pi*rand(1, 1);
    phi(5) = 2*pi*rand(1, 1);
    phi(6) = phi(4) + 2*phi(5);
    phis = phi';

    % create random signal
    syms t;
    x = cos(2*pi*lambda(1)*t + phi(1));
    for i = 2:6
        x = x + cos(2*pi*lambda(i)*t + phi(i));
    end
    x = matlabFunction(x);
    
    if display ~= 0
        figure();
        plot(1:1:N, x(0:1:N-1));
        title("Time series");
    end
    
    samples = x(0:1:N-1);
end