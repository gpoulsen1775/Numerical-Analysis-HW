
%HOMEWORK 7: Probelm 1


% Boundary Condititons  
a = 0; % u(0) = 0
b = 0; % u'(1) = 0

%ODE & Exact sol. 
f = @(x) sin(pi*x);
u_exact = @(x) -sin(pi*x)/(pi*pi);

%Num of Points, xVec, dt
N = 32;
x = linspace(0,1,N)';
h = x(2) - x(1);

% Right hand side
F = f(x(2:end-1));

% Fix up for boundary conditions
F(1) = F(1) - a/h^2;
F(end) = F(end) - (b*h + a)/h^2;

% Create matrix skipping first & last elems. 
A = zeros(N-2,N-2);

%Fill Corners
A(1,1) = -2/h^2;
A(1,2) = 1/h^2;
A(end,end) = -2/h^2;
A(end,end-1) = 1/h^2;

%Fill Rest of the Matrix
for ii = 2:(N-3)
    A(ii,ii) = -2/h^2;
    A(ii,ii-1) = 1/h^2;
    A(ii,ii+1) = 1/h^2;
end

% Solve system
UU = linsolve(A,F);
U = [a; UU; b];

%Plot 
figure(1); clf; hold on;
plot(x,U,'linewidth',4);
plot(x,u_exact(x),'g--', 'linewidth',4);