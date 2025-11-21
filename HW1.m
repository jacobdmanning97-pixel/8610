%Problem 2
% Initial conditions
x(1) = 11/2;
x(2) = 61/11;

% Number of iterations
N = 100;

% Recurrence relation
for k = 2:N
    x(k+1) = 111 - (1130 - 3000/x(k-1))/x(k);
end

% Exact solution
a = 0; b = 1; c = 1;
for k = 1:N
    exact(k) = (100^(k+1)*a + 6^(k+1)*b + 5^(k+1)*c) / (100^k*a + 6^k*b + 5^k*c);
end

% Plotting
figure
plot(1:N+1, x, 'r', 1:N, exact, 'b')
legend('Recurrence', 'Exact solution')
xlabel('k')
ylabel('x_k')
title('Comparison of Recurrence and Exact Solution')

%Problem 3
r = roots(wilkinson24_monomial_coeffs)

ls = [];
for j= 1 : 5
    for i = 1:5
        a = wilkinson24_monomial_coeffs(6+i);
        k = (19-j)^(18-i)*a/(factorial(18-j)*factorial(5+j));
        ls = [ls,abs(k)];
    end
end

max(ls)