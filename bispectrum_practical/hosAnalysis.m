function [powershift, bispPa, bispNo, frequency, bispDirect, waxis] = hosAnalysis(samples, M, maxLag, display)    
    % PSD
    powershift = powerspectrum(samples, M, 128, display);
    
    % bispectrum using HOSA toolbox
    NFFT = 2*maxLag + 1;
    [bispPa, frequency] = bispeci(samples, maxLag, M, 0, 'unbiased', NFFT, -1, display); % rect 
    [bispNo] = bispeci(samples, maxLag, M, 0, 'unbiased', NFFT, 0, display);  % parzen

    %[bispPa, bispNo, frequency] = bispAlternative(samples, maxLag, display);
    
    % bispectrum direct
    NFFT = M;
    [bispDirect, waxis] = bispecd(samples, NFFT, -1, M, 0, display);
end


function [bispPa, bispNo, frequency] = bispAlternative(samples, maxLag, display)
    % Bispectrum estimation using another toolbox "bisp3cum"
    
    [bispPa, frequency, cum3alt] = bisp3cum(samples, 1, maxLag, 'pa', 'unbiased', 0);
    if display ~= 0
        %figure();
        %mesh(-maxLag:maxLag, -maxLag:maxLag, cum3alt);
        %title("Third order cumulant bisp3cum");
        figure();
        subplot(211);
        hold on;
        plot(frequency(maxLag+1:end), frequency(maxLag+1:end), 'color', 'red');
        contour(frequency(maxLag+1:end), frequency(maxLag+1:end), ...
            abs(bispPa(maxLag+1:end, maxLag+1:end))), colorbar;
        title("Bispectrum Indirect Parzen");
        subplot(212);
        mesh(frequency(maxLag+1:end), frequency(maxLag+1:end), abs(bispPa(maxLag+1:end,maxLag+1:end)));
        colorbar;
    end
    
    [bispNo, frequency] = bisp3cum(samples, 1, maxLag, 'none', 'unbiased', 0);
    if display ~= 0
        figure();
        subplot(211);
        hold on;
        plot(frequency(maxLag+1:end), frequency(maxLag+1:end), 'color', 'red');
        contour(frequency(maxLag+1:end), frequency(maxLag+1:end), ...
            abs(bispNo(maxLag+1:end, maxLag+1:end))), colorbar;
        title("Bispectrum Indirect No window");
        subplot(212);
        mesh(frequency(maxLag+1:end), frequency(maxLag+1:end), abs(bispNo(maxLag+1:end, maxLag+1:end)));
    end
end