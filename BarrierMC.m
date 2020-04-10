%% Parameters
S0 = 100;
K = 100;
B = 130;
q = 0.05;
r = 0.03;
T = 0.5;
alpha = 0.35;

% N and M
N = 1000000;
M = 500;

% Specify time paritions between terminal time points 
dt = T/M;
tm = 0:dt:T;
m = length(tm)-1;
dtm = diff(tm);

% Create vectors to store the original and the antithetic sample
Sample = zeros(m+1,1);
S_antithetic = zeros(m+1,1);

% Create vector to store the payoff for both samples
payoff = zeros(N,2);

% Initial price of sample points 
Sample(1) = 100;
S_antithetic(1) = 100;


% Create a for loop for discretizing the Stochastic Process
for i=1:N
    
    
    for k=1:m
        
        z = randn;
        
        % Discretizing original sample 
        Sample(k+1) = Sample(k) + (r-q)*Sample(k)*dtm(k) + 0.25*exp(-tm(k))*(100^alpha)*(Sample(k)^(1-alpha))*...
            sqrt(dtm(k))*z;
        
        % Discretizing antithetic sample
        S_antithetic(k+1) = S_antithetic(k) + (r-q)*S_antithetic(k)*dtm(k) - 0.25*exp(-tm(k))*(100^alpha)*(S_antithetic(k)^(1-alpha))*...
            sqrt(dtm(k))*z;
        
    end
    
  % Creating a for-loop to incorporate the barrier Ceiling for original
  % sample
    if any(Sample >= 130)
        payoff(i,1) = 0;
    else
        payoff(i,1) = exp(-T*r)*max(Sample(m+1)-K,0);
    end
       
    % Creating a for-loop to incorporate the barrier Ceiling for antithetic
  % sample
    if any(S_antithetic >= 130)
        payoff(i,2) = 0;
    else
        payoff(i,2) = exp(-T*r)*max(S_antithetic(m+1)-K,0);
    end
    
    
    
end

% Compute the Barrier Option price
Barrier_mc = mean(mean(payoff));

% Compute the standard error of the sample
SE_C = sqrt((sum(mean(payoff,2).^2)-N*Barrier_mc^2)/(N*(N-1)));

% Compute the confidence interval.
CI_low = Barrier_mc - 1.96*SE_C;
CI_up = Barrier_mc + 1.96*SE_C;

% Plotting samples

figure(1)
plot(tm,Sample)
xlabel('Time')
ylabel('Stock Price')
title('Stock Price path in Monte Carlo using Euler Discretisation')
legend('Original Sample Stock Price')

figure(2)
plot(tm,S_antithetic)
xlabel('Time')
ylabel('Stock Price')
title('Stock Price path in Monte Carlo using Euler Discretisation')
legend('Antithetic Sample Stock Price')


