% Function to compute the product between a tri-diagonal matrix (d u l) and a
% vector x
function y = tridiag_prod(d,u,l,x)
N = length(d);
y = zeros(N,1);
for j=1:N
        if(j==1)
            y(j) = d(j)*x(j) + u(j)*x(j+1);
        elseif(j<N)
            y(j) = l(j-1)*x(j-1) + d(j)*x(j) + u(j)*x(j+1);
        else
            y(j) = l(j-1)*x(j-1) + d(j)*x(j);
        end
end
end