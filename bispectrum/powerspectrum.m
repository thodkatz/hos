function [psd] = powerspectrum(samples, M, maxLag, display)
    if M == maxLag
        M = M + 1; 
    end    

    cum2 = cumest(samples, 2, maxLag, M, 0, 'unbiased');
    
    if display ~= 0
        figure();
        plot(0:maxLag, cum2(maxLag+1:end)); % symmetric
        title("Autocorrelation cumest");
    end

    % power spectrum
    fs = 1;
    fshift = (-maxLag:maxLag)*(fs/(2*maxLag+1));
    psd = abs(fftshift(fft(cum2)));
    
    if display ~= 0
        figure();
        plot(fshift, psd);
        grid on;
        title("PSD");
        xlabel("frequency (Hz)");
        DragDataTip
    end
end