clear;
clc;
close all;

display = 0;
N = 8192;
samples = genSamples(N, display);

maxLag = 128;
M = 256;

% repeat for different M = 512 and M = 128
hosAnalysis(samples, 128, maxLag, 0);
hosAnalysis(samples, 256, maxLag, 0);
hosAnalysis(samples, 512, maxLag, 0);

% repeat for 50 realizations and get the mean
numIter = 10;
psd = zeros(numIter, 2*maxLag + 1);
bispPa = zeros(numIter, 2*maxLag + 1, 2*maxLag + 1);
bispNo = zeros(numIter, 2*maxLag + 1, 2*maxLag + 1);
bispDirect = zeros(numIter, M, M);
for k = 1:numIter
    samples = genSamples(N, display);
    [psd(k,:),bispPa(k,:,:), bispNo(k,:,:), frequency, bispDirect(k,:,:), waxis] = ...
        hosAnalysis(samples, M, maxLag, 0);
end

bispNoVariance = zeros(2*maxLag+1, 2*maxLag+1);
bispPaVariance = zeros(2*maxLag+1, 2*maxLag+1);
bispDirectVariance = zeros(M, M);
for i = 1:2*maxLag+1
   for j = 1:2*maxLag+1
        bispPaVariance(i,j) = var(bispPa(:,i,j));
        bispNoVariance(i,j) = var(bispNo(:,i,j));
   end
end

for i = 1:M
    for j = 1:M 
        bispDirectVariance(i,j) = var(bispDirect(:,i,j));
    end
end

psdMean = mean(psd, 1);
bispPaMean = mean(bispPa, 1);
bispNoMean = mean(bispNo, 1);
bispDirectMean = mean(bispDirect, 1);
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
    title("PSD variance");
end

if display ~= 0
    figure();
    subplot(211);
    contour(frequency, frequency, abs(bispPaMean));
    title("Bispectrum bisp3cum parzen");
    subplot(212);
    mesh(frequency(maxLag+1:end), frequency(maxLag+1:end), ...
        abs(bispPaMean(maxLag+1:end,maxLag+1:end)));
        
    figure(); 
    mesh(frequency(maxLag+1:end), frequency(maxLag+1:end), ...
        abs(bispPaVariance(maxLag+1:end,maxLag+1:end)));
    title("Bispectrum Parzen variance");
    
    figure();
    subplot(211);
    contour(frequency, frequency, abs(bispNoMean));
    title("Bispectrum bisp3cum no window");
    subplot(212);
    mesh(frequency(maxLag+1:end), frequency(maxLag+1:end), ...
        abs(bispNoMean(maxLag+1:end, maxLag+1:end)));
    
    figure(); 
    mesh(frequency(maxLag+1:end), frequency(maxLag+1:end), ...
        abs(bispNoVariance(maxLag+1:end,maxLag+1:end)));
    title("Bispectrum No window variance");
    
    figure();
    subplot(211);
    contour(waxis, waxis, abs(bispDirectMean));
    title("Bispectrum bisp3cum Direct");
    subplot(212);
    mesh(waxis(M/2:end), waxis(M/2:end), ...
        abs(bispDirectMean(M/2:end, M/2:end)));
    
    figure(); 
    mesh(waxis(M/2:end), waxis(M/2:end), ...
        abs(bispDirectVariance(M/2:end,M/2:end)));
    title("Bispectrum Direct variance");
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
    phi(6) = phi(4) + phi(5);
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