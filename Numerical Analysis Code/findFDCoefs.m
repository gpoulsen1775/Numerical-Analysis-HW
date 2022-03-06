function coefs = findFDCoefs(k, xbar, xvals)

% Logical Check
N = length(xvals);
if (N < k)
    error("Not enough x values!");
end

% Make Matrix
mat = zeros(N,N);

% i will loop through each row
for i = 1:N
    
    % Develop Taylor Seires Leading Coefficient 
    fac = 1.0/factorial(i - 1);
    
    % j will loop through each column
    for j = 1:N
        mat(i,j) = fac * (xvals(j) - xbar)^(i - 1);
    end
end

% RHS has zeros everywhere except for the derivative.
% Make Result Vector: k+1 cuz we start from 0th derivative
rhs = zeros(N,1);
rhs(k+1) = 1.0;

% Solve the system.
coefs = linsolve(mat, rhs)';

end