function [powershift, bispPa, bispNo, frequency] = hosAnalysis(samples, M, maxLag, display)    
    % autocorrelation hosa toolbox
    cum2 = cumest(samples, 2, maxLag, 0, 0, 'biased');
    % cumest dont work well with samples seg and unbiased

    % normalizing with the variance
    % should the autocorrelation be normalized for power spectrum?
    %cum2 = cum2 ./ var(samples);

    if display ~= 0
        figure();
        %plot(-maxLag:maxLag, cum2);
        plot(0:maxLag, cum2(maxLag+1:end)); % symmetric
        title("Autocorrelation cumest");
    end

    % autocorrelation via xcorr
    %acf = xcorr(samples, maxLag, 'unbiased');
    if display ~= 0
        %figure();
        %plot(-maxLag:maxLag, acf);
        %plot(0:maxLag, acf(maxLag+1:end));
        %title("Autocorrelation xcorr");
    end

    % power spectrum
    fs = 1;
    fshift = (-maxLag:maxLag)*(fs/(2*maxLag+1)); % zero-centered frequency range
    
    %powershift = abs(fftshift(fft(acf)));
    if display ~= 0
        %figure();
        %plot(fshift,powershift);
        %plot(fshift(maxLag+1:end), powershift(maxLag+1:end));
        %title("Power spectral density autocorr");
    end

    powershift = abs(fftshift(fft(cum2)));
    
    if display ~= 0
        figure();
        plot(fshift, powershift);
        title("Power spectral density cumest");
    end
    
    % test for fft. Does it matter getting only the half positive vs all
    %powershift = abs(fftshift(fft(cum2(maxLag:2*maxLag-1))));
    
    if display ~= 0
        %figure();
        %plot(-maxLag/2:maxLag/2-1, powershift);
        %title("Power spectral density cumest 2");
    end
    
    % specs
    %maxLag = 64; % homework
    
    % third order cumulant
    %for k = -maxLag:maxLag
    %cum3(:,k+maxLag+1) = cumest(samples, 3, maxLag, 0, 0, 'unbiased', k); % won't work for M as sample length
    %end

    % should the autocorrelation be normalized with kyrtosis?
    %kyrtosis = cum3(maxLag+1,maxLag+1) % c3(0,0)
    %cum3 = cum3 ./ kyrtosis;

    if display ~= 0
        %figure();
        %subplot(211);
        %mesh(-maxLag:maxLag, -maxLag:maxLag, cum3);
        %title("Third order cumulant cumest");
        %subplot(212);
        %contour(-maxLag:maxLag, -maxLag:maxLag, cum3, 8);
    end
    
    % bispectrum indirect
    NFFT = 2*maxLag + 1;
    %bispeci(samples, maxLag, M, 0, 'unbiased', NFFT, -1, display); % rect 
    %bispeci(samples, maxLag, M, 0, 'unbiased', NFFT, 0, display);  % parzen

    % bispectrum alternative
    [bispPa, frequency, cum3alt] = bisp3cum(samples, 1, maxLag, 'pa', 'unbiased', display);
    if display ~= 0
        figure();
        mesh(-maxLag:maxLag, -maxLag:maxLag, cum3alt);
        title("Third order cumulant bisp3cum");
        figure();
        subplot(211);
        contour(frequency, frequency, abs(bispPa));
        title("Bispectrum bisp3cum parzen");
        subplot(212);
        mesh(frequency(maxLag+1:end), frequency(maxLag+1:end), abs(bispPa(maxLag+1:end,maxLag+1:end)));
    end
    
    [bispNo, frequency] = bisp3cum(samples, 1, maxLag, 'none', 'unbiased', display);
    if display ~= 0
        figure();
        subplot(211);
        contour(frequency, frequency, abs(bispNo));
        title("Bispectrum bisp3cum no window");
        subplot(212);
        mesh(frequency(maxLag+1:end), frequency(maxLag+1:end), abs(bispNo(maxLag+1:end, maxLag+1:end)));
    end
    
    % find the index of the array with the peak?
    % bispectrum magnitude give us the frequencies that are coupled. How
    % they are coupled, the outcome (addition?) can be shown via the phase
    % of the bispectrum?
    
    % bispectrum direct
    NFFT = M;
    bispecd(samples, NFFT, -1, M, 0, display);
end