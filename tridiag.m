function x = tridiag(d,u,l,b)
n = length(b);
x = zeros(n,1);
for j=2:n
    d(j) = d(j) - l(j)*u(j-1)/d(j-1);
    b(j) = b(j) - l(j)*b(j-1)/d(j-1);
end
x(n) = b(n)/d(n);
for j=n-1:-1:1
    x(j) = (b(j) - u(j)*x(j+1))/d(j);
end
end