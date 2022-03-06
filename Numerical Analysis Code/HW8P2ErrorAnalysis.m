
%%%%%%%%%%%%%%%%%%%Set Up%%%%%%%%%%%%%%%%%%%%%%

f = @(x,y) (0 * x).* (0*y);
u_exact = @(x,y) x + y;

Narray = [32, 64, 128, 256];
Harray = [1/32, 1/64, 1/128, 1/256];

Error = 0*Narray;

for i = 1:4 

N = Narray(1,i);
x = linspace(0,1,N)';
y = linspace(0,1,N)';
m = N - 2; 
h = x(2) - x(1);
[xx,yy] = meshgrid(x,y);

% Find F
F = f(xx(2:end-1,2:end-1),yy(2:end-1,2:end-1));
F = F(:); % This turns F into a column vector that is row-wise ordering

%%%%%%%%%%%%Boundary Conditions%%%%%%%%%%%%%%%%%

% For Row One
F(1:m) = F(1:m) - u_exact(x(2:end-1),y(1)) / h^2;

% For Row m
F((end-m+1):end) = F((end-m+1):end) - u_exact(x(2:end-1),y(end)) / h^2;

% Every 1st point on each row
F(1:m:end) = F(1:m:end) - u_exact(x(1),y(2:end-1)) / h^2;

% Every Last point on each row 
F(m:m:end) = F(m:m:end) - u_exact(x(end),y(2:end-1)) / h^2;

%%%%%%%%%%%%%%%%Create Matrix A%%%%%%%%%%%%%%%%%

% Create matrix
% Make it sparse
I = eye(m);
e = ones(m,1);
T = spdiags([4*e -20*e 4*e], [-1 0 1], m, m);
R = spdiags([e 4*e e],[-1 0 1], m, m);
S = spdiags([e e], [-1 1], m, m);
A = (kron(I, T) + kron(S, R))/h^2;

UU = A\F;

%%%%%%%%%%%%%%%%%Error Analysis%%%%%%%%%%%%%%%%%%%%%%%

exactSol = u_exact(xx(2:end-1,2:end-1),yy(2:end-1,2:end-1));
a = reshape(exactSol, m*m, 1);
E = abs(UU - a);
tao = max(A * E);

Error(1,i) = tao;

end

%%%%%%%%%%%%%%%Error Plotting%%%%%%%%%%%%%%%%%%%%%%%%
plot(Harray,Error);
hold on;
fplot(@(x) x^4, [0,.03]);