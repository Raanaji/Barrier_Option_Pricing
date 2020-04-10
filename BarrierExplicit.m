% Explicit Scheme to solve Barrier Call Option
%For further instructions, please refer the CODE DOCUMENTATION
 
%% The Parameters 
% Strike Prices

K = 100;

S0 = 100;

% Set the minimal and maximal stock prices

Smin = 0;

Smax = 130;

% The number of points in stock direction

N = 1000;

S1 = linspace(Smin,Smax,N+1)';

% The length of the stock price interval

dS = S1(2) - S1(1);

% S stores all the prices except boundary points

S = S1(2:N);

 

% Maturity of option

T = 0.5;

% The number of points in time dimension

M = 50000;

tau = linspace(0,T,M+1);

% The length of the time interval

dtau = tau(2) - tau(1);

 

% Interest rate

r = 0.03;

% Dividend yield

q = 0.05;

a = 0.35; % Volatility Alpha

sigma = zeros(N-1,M);

for k=1:M

    for j=1:N-1

sigma(j,k)=0.25*exp(-tau(k))*((100/S(j))^a);

    end

end

% alpha and beta in the discretization scheme

alpha = 0.5*(sigma.^2).*(S.^2)*dtau/(dS^2);

beta = (r - q)*S*dtau/(2*dS);

 

% lower, main and upper diagonal of the tridiagonal matrix for the explicit

% finite difference scheme

l = alpha - beta;

d = 1 - r*dtau - 2*alpha;

u = alpha + beta;

 


% Vector to store the option prices at time k

Vold = max(S - K,0);

% Vector to store the option prices at time k+1

Vnew = Vold;

for k=1:M
    % Boundary condition for Barrier Call option
    
   
    boundary = [l(1)*0;zeros(N-3,1);u(N-1)*0];
    
    % Explicit iteration scheme
    for j=1:N-1
        if(j==1)
            Vnew(j) = d(j,k)*Vold(j) + u(j,k)*Vold(j+1);
        elseif(j<N-1)
            Vnew(j) = l(j,k)*Vold(j-1) + d(j,k)*Vold(j) + u(j,k)*Vold(j+1);
        else
            Vnew(j) = l(j,k)*Vold(j-1) + d(j,k)*Vold(j);
        end
    end
    % Update the vectors from time k to time k+1
    Vold = Vnew + boundary;
end

 

% Interpolation to find the call price when S0 = 100

call_fdm = interp1(S,Vold,S0);
