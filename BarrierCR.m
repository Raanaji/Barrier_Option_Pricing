% Crank Nicholson Scheme to solve Barrier Call Option
%For further instructions, please refer the CODE DOCUMENTATION
 
%% The Parameters 
% Strike Price

K = 100;

S0 = 100;

% Set the minimal and maximal stock prices

Smin = 0;

Smax = 130;

% The number of points in stock direction

N = 100;

S1 = linspace(Smin,Smax,N+1)';

% The length of the stock price interval

dS = S1(2) - S1(1);

% S stores all the prices except boundary points

S = S1(2:N);

 

% Maturity of option

T = 0.5;

% The number of points in time dimension

M = 5000;

tau = linspace(0,T,M+1);

% The length of the time interval

dtau = tau(2) - tau(1);

 

% Interest rate

r = 0.03;

% Dividend yield

q = 0.05;

a = 0.35; % Volatility Alpha

% Discretized Volatility

sigma = zeros(N-1,M);

for k=1:M

    for j=1:N-1
        %-T+tau(k) 

sigma(j,k)=0.25*exp(-tau(k))*((100/S(j))^a);

    end

end

% alpha and beta in the discretization scheme

alpha = 0.5*(sigma.^2).*(S.^2)*dtau/(dS^2);

beta = (r - q)*S*dtau/(2*dS);

 

% lower, main and upper diagonal of the tridiagonal matrix for the implicit

% finite difference scheme

l = -alpha + beta;

d = 1 + r*dtau + 2*alpha;

u = -(alpha + beta);

 

% lower, main and upper diagonal of the tridiagonal matrix for the explicit

% finite difference scheme

lE = alpha - beta;

dE = 1 - r*dtau - 2*alpha;

uE = alpha + beta;


% Vector to store the option prices at time k

% Option payoff at maturity

Vold = max(S - K,0);

% Vector to store the option prices at time k+1

Vnew = Vold;

temp = Vold;

for k=1:M
    % Boundary condition for Barrier Call option at time k+1
    
    boundary_I = [-l(1)*0;zeros(N-3,1);-u(N-1)*0];
    
    % Boundary condition for Barrier Call option at time k
   
    boundary = [l(1)*0;zeros(N-3,1);uE(N-1)*0]; 
    
    % To make sure 2nd order convergence in time, a couple of implicit
    % scheme is implemented before switching to crank nicolson scheme
    
    if(k<3)
        Vnew = tridiag(d(:,k),u(:,k),l(:,k),Vold+boundary_I);
    else    
        
    % Crank Nicolson iteration scheme
    % The correct scheme is that A+I, hence instead of dE, we need
    % to use dE+1
    % Explicit iteration scheme
    for j=1:N-1
        if(j==1)
            temp(j) = (dE(j,k)+1)*Vold(j) + uE(j,k)*Vold(j+1);
        elseif(j<N-1)
            temp(j) = lE(j,k)*Vold(j-1) + (dE(j,k)+1)*Vold(j) + uE(j,k)*Vold(j+1);
        else
            temp(j) = lE(j,k)*Vold(j-1) + (dE(j,k)+1)*Vold(j);
        end
    end
 
    Vnew = tridiag(d(:,k)+1,u(:,k),l(:,k),temp+boundary);
    end
    % Update the vectors from time k to time k+1
    Vold = Vnew; 
end

% Interpolation to find the call price when S0 = 100

call_fdm = interp1(S,Vold,S0);