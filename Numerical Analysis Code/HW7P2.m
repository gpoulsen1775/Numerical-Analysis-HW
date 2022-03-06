
%%Original Code Modified To Plot Sol. W/ 4th Order Methods

a = 0; % u(0) = 0
b = 0; % u(1) = 0

% ODE & Sol'n
f = @(x) sin(pi*x);
u_exact = @(x) -sin(pi*x)/(pi*pi);

% Num of pts, xVec, dt
n = 32;
x = linspace(0,1,n)';
h = x(2) - x(1);

% Right hand side
F = f(x(2:end-1));

% Specify F for Row 1 & 2
F(1) = F(1) - a/(12*h^2);
F(2) = F(2) + a/(12*h^2);

% Specify F for Row end-1
F(end-1) = F(end) + b/(12*h^2);

% Create matrix skipping first & last elems. 
A = zeros(n-2,n-2);

%Set Row 1 W/ F. Diff. Stencil
A(1,1) = (1/(12*h^2))*-15;
A(1,2) = (1/(12*h^2))*-4;
A(1,3) = (1/(12*h^2))*14;
A(1,4) = (1/(12*h^2))*-6;
A(1,5) = (1/(12*h^2))*1;

%Set Row 2
A(2,1) = (1/(12*h^2))*16;
A(2,2) = (1/(12*h^2))*-30;
A(2,3) = (1/(12*h^2))*16;
A(2,4) = (1/(12*h^2))*-1;

%Set Row end-1 
A(end-1,end) = (1/(12*h^2))*16;
A(end-1,end-1) = (1/(12*h^2))*-30;
A(end-1,end-2) = (1/(12*h^2))*16;
A(end-1,end-3) = (1/(12*h^2))*-1;

%Set Row end W/ Backwards Diff. Stencil
A(end,end) = (1/(12*h^2))*24;
A(end,end-1) = (1/(12*h^2))*-60;
A(end,end-2) = (1/(12*h^2))*48;
A(end,end-3) = (1/(12*h^2))*-12;

%Fill Rest of Matrix W/ Center Diff. W/ all 5 pts. 
for ii = 3:(n-4)
    
    A(ii,ii) = (1/(12*h^2))*-30;
    A(ii,ii+1) = (1/(12*h^2))*16;
    A(ii,ii-1) = (1/(12*h^2))*16;
    A(ii,ii+2) = (1/(12*h^2))*-1;
    A(ii,ii-2) = (1/(12*h^2))*-1;
end

% Solve system
UU = linsolve(A,F);
U = [a; UU; b];

%Plot Sol'n Against Apporx. 
figure(1); clf; hold on;
plot(x,U,'b','linewidth',4);
plot(x,u_exact(x),'g--', 'linewidth',4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Now We Check the Error of Our Methods:

%Number of Points used & dt used
N = [32, 64, 128, 256];
H = [1/32; 1/64; 1/128; 1/256];

%Global Error Vector
Err1 = zeros(32,1);
Err2 = zeros(64,1);
Err3 = zeros(128,1);
Err4 = zeros(256,1);

for i = 1:4 
    
a = 0; % u(0) = 0
b = 0; % u(1) = 0

% ODE & Sol'n
f = @(x) sin(pi*x);
u_exact = @(x) -sin(pi*x)/(pi*pi);

% Num of pts, xVec, dt
n = N(1,i);
x = linspace(0,1,n)';
h = x(2) - x(1);

% Right hand side
F = f(x(2:end-1));

% Specify F for Row 1 & 2
F(1) = F(1) - a/(12*h^2);
F(2) = F(2) + a/(12*h^2);

% Specify F for Row end-1
F(end-1) = F(end) + b/(12*h^2);

% Create matrix skipping first & last elems. 
A = zeros(n-2,n-2);

%Set Row 1 W/ F. Diff. Stencil
A(1,1) = (1/(12*h^2))*-15;
A(1,2) = (1/(12*h^2))*-4;
A(1,3) = (1/(12*h^2))*14;
A(1,4) = (1/(12*h^2))*-6;
A(1,5) = (1/(12*h^2))*1;

%Set Row 2
A(2,1) = (1/(12*h^2))*16;
A(2,2) = (1/(12*h^2))*-30;
A(2,3) = (1/(12*h^2))*16;
A(2,4) = (1/(12*h^2))*-1;

%Set Row end-1 
A(end-1,end) = (1/(12*h^2))*16;
A(end-1,end-1) = (1/(12*h^2))*-30;
A(end-1,end-2) = (1/(12*h^2))*16;
A(end-1,end-3) = (1/(12*h^2))*-1;

%Set Row end W/ Backwards Diff. Stencil
A(end,end) = (1/(12*h^2))*24;
A(end,end-1) = (1/(12*h^2))*-60;
A(end,end-2) = (1/(12*h^2))*48;
A(end,end-3) = (1/(12*h^2))*-12;

%Fill Rest of Matrix W/ Center Diff. W/ all 5 pts. 
for ii = 3:(n-4)
    
    A(ii,ii) = (1/(12*h^2))*-30;
    A(ii,ii+1) = (1/(12*h^2))*16;
    A(ii,ii-1) = (1/(12*h^2))*16;
    A(ii,ii+2) = (1/(12*h^2))*-1;
    A(ii,ii-2) = (1/(12*h^2))*-1;
end

% Solve system
UU = linsolve(A,F);
U = [a; UU; b];
    
%Accumulated Error of Sol. 
err = zeros(n,1);

%Find Accumulated Error for choice of n
for j = 1: n
    
    err(j,1) = abs(U(j)-f((j-1)*h));
end

if (i == 1)
    Err1 = U.*err;
end

if (i == 2)
    Err2 = U.*err;
end

if (i == 3)
    Err3 = U.*err;
end

if (i == 4)
    Err4 = U.*err;
end

end

x1 = linspace(0,1,32)';
x2 = linspace(0,1,64)';
x3 = linspace(0,1,128)';
x4 = linspace(0,1,256)';

% Each respective Err is A.*E = -tao
% So taking the abs(-tao) yields truncation 
% error over solution. So we plot this for 
% error examination. 

figure(2);
hold on;
scatter(x1,abs(Err1));
fplot(@(x) x^3, [0,.5]);

figure(3); hold on;
scatter(x2,abs(Err2));
fplot(@(x) x^3, [0,.5]);

figure(4); hold on;
scatter(x3,abs(Err3));
fplot(@(x) x^3, [0,.5]);

figure(5); hold on;
scatter(x4,abs(Err4));
fplot(@(x) x^3, [0,.5]);

