function decay1(tol)

% decay1.m MODIFIED
% Sample code for solving a system of ODEs in matlab.
% Solves  system arising from decay process  A --> B --> C --> D.
%
% Adopted from  http://www.amath.washington.edu/~rjl/fdmbook/chapter8  (2007)

% Declare some global values. These parameters are used to calculate total
% number of function evaluations and are constants used in evaluating the
% ODE.
global fcnevals K1 K2 K3

% decay rates:
K1 = 1;
K2 = 2;
K3 = 3; 

t0 = 0;                      % initial time
u0 = [1 0 0 0];                % initial data for u(t) as a vector
tfinal = 4;                  % final time
fcnevals = 0;                % counter for number of function evaluations

% solve ode:
options = odeset('AbsTol',tol,'RelTol',tol);
odesolution = ode113(@f,[t0 tfinal],u0,options);

% Plot the solution
clf
hold on
t = linspace(0, tfinal, 500);
u = deval(odesolution, t);
plot(t,u(1,:),'b') % Component 1
plot(t,u(2,:),'k') % Component 2
plot(t,u(3,:),'r') % Component 3
plot(t,u(4,:),'g') % Component 3
legend('u1','u2','u3','u4')
axis([0 tfinal -.1 1.1])
title('u(t) as a function of time')


%----------------------------------

function f = f(t,u)
% Grab the global variables from above
global fcnevals K1 K2 K3

% Evaluate the function
f1 = -K1*u(1);
f2 = K3*u(4) - K2*u(2); 
f3 = K2*u(2);
f4 = K1*u(1) -K3*u(4);
f = [f1; f2; f3; f4]; % Form final vector of function evaluations.

fcnevals = fcnevals + 1;

