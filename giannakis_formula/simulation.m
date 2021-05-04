clc;
clear;
close;

% input
N = 2048;
v = exprnd(1, [1 N]);

% output MA(5)
q = 5;
coeff = [1.0 0.93 0.85 0.72 0.59 -0.10];
x = conv(v, coeff, 'same');

% calculate skewness
skewness = sum((v - mean(v)).^3)/(N-1)*std(v)^3;

% 3rd order cumulant
maxLag = 20;
M = 64;
nsamp = 2*maxLag + 1;
cum3 = zeros(nsamp, nsamp);
for k = -maxLag:maxLag
  cum3(:,k+maxLag+1) = cumest(x, 3, maxLag, M, 0, 'biased', k);
end
mesh(-maxLag:maxLag, -maxLag:maxLag, cum3)
title("Third order cumulant of output")

% Estimate impulse response with Giannakis formula
order = [q q-2 q+3]; % include sub and sup estimation
h_hat = zeros(N,3);
for i = 1:3
h_hat(:, i) = [cum3(order(i),maxLag+1:maxLag+1+order(i)) ./ ...
    cum3(order(i),maxLag+1) zeros(1, N-order(i)-1)];
end

normal = h_hat(1:order(1) + 1,1)
sub    = h_hat(1:order(2) + 1,2)
sup    = h_hat(1:order(3) + 1,3)

% Estimate system output with h_hat
