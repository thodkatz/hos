clc;
clear;
close all;

% samples
N = 2^11
snr = -5: 5: 30;
numIters = 10;
nrms = zeros(numel(snr), numIters);

for i = 1:numIters
    fprintf("%d\n", i)
    nrms(:,i) = per_realization(N,snr);
end

figure
plot(snr, mean(nrms, 2)')
xlabel("SNR (dB)")
ylabel("Average NRMSE")

function nrms = per_realization(N, snr)
    % input
    v = exprnd(1, [1 N]);
    v = v - mean(v);

    % calculate skewness
    skewness = sum((v - mean(v)).^3)/(N-1)*std(v)^3;

    % output MA(5)
    q = 5;
    h = [1.0 0.93 0.85 0.72 0.59 -0.10];
    x = conv(v, h, 'same');
    %x = x - mean(x);

    % output estimated for three different order and for three different SNRs
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
    
    % get the rnms for the q order
    nrms = nrms(1,:);
    
    display_snr = 0;
    if display_snr ~= 0
        figure
        plot(snr, nrms)
        xlabel("SNR (dB)")
        ylabel("NRMSE")
    end
end
