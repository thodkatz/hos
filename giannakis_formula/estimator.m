function [x_est,rmse,nrms] = estimator(v,x,q,N, display)
    % 3rd order cumulant
    maxLag = 20;
    M = 64;
    nsamp = 2*maxLag + 1;
    cum3 = zeros(nsamp, nsamp);
    for k = -maxLag:maxLag
      cum3(:,k+maxLag+1) = cumest(x, 3, maxLag, M, 0, 'unbiased', k);
    end
    
    if display ~= 0
        figure
        mesh(-maxLag:maxLag, -maxLag:maxLag, cum3)
        title("Third order cumulant of output")
    end
    
    % Estimate impulse response with Giannakis formula
    order = [q q-2 q+3]; % include sub and sup estimation
    h_hat = zeros(N,3);
    for i = 1:3
    h_hat(:, i) = [cum3(order(i),maxLag+1:maxLag+1+order(i)) ./ ...
        cum3(order(i),maxLag+1) zeros(1, N-order(i)-1)];
    end

    normal = h_hat(1:order(1) + 1,1);
    sub    = h_hat(1:order(2) + 1,2);
    sup    = h_hat(1:order(3) + 1,3);

    % Estimate system output with h_hat
    x_est = zeros(N,3);
    for i = 1:3
        x_est(:,i) = conv(v, h_hat(1:order(i) + 1,1), 'same');
        %x_est( = x_est - mean(x_est);
    end
    
    if display ~= 0
        figure
        hold on
        plot(x_est(:,1), 'color', 'red')
        plot(x, 'color', 'blue')
        title("Order " + order(1))
    end
    
    if display ~= 0
        figure
        hold on
        plot(x_est(:,2), 'color', 'red')
        plot(x, 'color', 'blue')
        title("Order " + order(2))
    end

    if display ~= 0
        figure
        hold on
        plot(x_est(:,3), 'color', 'red')
        plot(x, 'color', 'blue')
        title("Order " + order(3))
    end
    
    % quantification errors
    rmse = zeros(1, 3);
    nrms = zeros(1, 3);
    for i = 1:3
        rmse(i) = sqrt(sum((x_est(:,i)' - x).^2)./N);
        nrms(i) = rmse(i) / (max(x) - min(x));
    end
end