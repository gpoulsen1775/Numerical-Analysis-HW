function [graph] = hw4p5()

% hw4p5: No Input; Ouput: Graph of Error of Leapfrog/Forward E. over time
% against a graph of Error of a linear order numerical Method (Forward E.
% Alone) over time. 

% The function evaluates a solution to the ODE u' = 2*u in two different ways
% (Leapfrog/Forward E. & Forward E.) & then compares the error of both
% methods using a graph in order to demonstrate that Leapfrog/Forward E.
% does not get reduced to a linear order approximation by implementing
% Forward E. for step 1. 

%////////////////////////////////////////////////////////////////////////

% First We Preform Forward Euler's & Leapfrog

% ODE, Parameters, Solution
lambda = 2;
dudt = @(u) lambda*u;
u0 = 2;
u = @(t) u0*exp(lambda*t);
delta = .001;

% Number of Time Steps
n = 10/delta;

% Create Solution Vector
sol1 = zeros(n,2);

% Initial Conditions
sol1(1,1) = 0; 
sol1(1,2) = u0;

% Preform Forward Euler's for U1
sol1(2,1) = sol1(1,1) + delta;
sol1(2,2) = sol1(1,2) + delta*dudt(sol1(1,2));

% Perform LeapFrog On the Rest
for i = 3: n

    sol1(i,1) = sol1(i-1,1) + delta; 
    sol1(i,2) = sol1(i-2,2) + 2*delta*dudt(sol1(i-1,2));

end

% We Now Calculate our Error
err1 = zeros(n,2); 

for i = 1: n

    err1(i,1) = sol1(i,1);
    err1(i,2) = abs( u(sol1(i,1)) - sol1(i,2));
end

% Plot Error of LeapFrog/Forward
figure(1);
plot(err1(:,1),err1(:,2),'g.');
hold on;

%////////////////////////////////////////////////////////////////////////

% Now we do the same for Forward Euler's alone
sol2 = zeros(n,2); 

% Initial Conditions
sol2(1,1) = 0;
sol2(1,2) = u0;

%Preform Forward Eulers
for i = 2: n
    
    sol2(i,1) = sol2(i-1,1) + delta;
    sol2(i,2) = sol2(i-1,2) + delta*dudt(sol2(i-1,2));
end

% Compute Error
err2 = zeros(n,2);

for i = 1: n
    
    err2(i,1) = sol2(i,1);
    err2(i,2) = abs( u(sol2(i,1)) - sol2(i,2));
end

% Plot Error For Forward Eulers Alone
figure(1);
plot(err2(:,1),err2(:,2),'r.');

title("Error of O(h) Method vs. Leapfrog/Euler's Error On [0,10]");

graph = figure(1);
end

