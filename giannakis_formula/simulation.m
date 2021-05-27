clc;
clear;
close all;

rng(1)

% samples
N = 2^11
q = 5;

[v,x] = create_output(N);
[x,y, nrms] = estimator(v,x,q,N,1);

nrms(1)
nrms(2)
nrms(3)

snr = -5:5:30;
per_realization_snr(N, q, snr);

numIters = 50;
%multiple_realizations_snr(N, q, numIters);

function multiple_realizations_snr(N,q,numIters)
    % snr multiple iterations
    snr = -5: 5: 30;
    nrms = zeros(numel(snr), numIters);

    for i = 1:numIters
        fprintf("%d\n", i)
        nrms(:,i) = per_realization_snr(N,q,snr);
    end

    figure
    plot(snr, mean(nrms, 2)')
    xlabel("SNR (dB)")
    ylabel("Average NRMSE")
    title(numIters + " realizations")

    figure
    plot(snr, std(nrms,0,2'))
    xlabel("SNR (dB)")
    ylabel("Standard deviation NRMSE")
    title(numIters + " realizations")
end

function nrms = per_realization_snr(N, q, snr)
    [v, x] = create_output(N);

    % output estimated for three different order and different SNRs
    x_est = zeros(N, 3, numel(snr));
    rmse = zeros(3, numel(snr));
    nrms = zeros(3, numel(snr));

    % loop for different snr values
    display_est = 0;
    for i = 1:numel(snr)
        % add gaussian noise
        noise = awgn(x, snr(i), 'measured');
        [x_est(:,:,i), rmse(:,i), nrms(:,i)] = estimator(v, noise, q, N, display_est);
    end
    
    % get the rnms for the q=5 order
    nrms = nrms(1,:);
    
    display_snr = 0;
    if display_snr ~= 0
        figure
        plot(snr, nrms)
        xlabel("SNR (dB)")
        ylabel("NRMSE")
    end
end

function [v, x] = create_output(N)
    % input
    v = exprnd(1, [1 N]);
    v = v - mean(v);

    % calculate skewness
    skewness = sum((v - mean(v)).^3)/(N-1)*std(v)^3;

    % output MA(5)
    h = [1.0 0.93 0.85 0.72 0.59 -0.10];
    x = conv(v, h, 'same');
    %x = x - mean(x);
end