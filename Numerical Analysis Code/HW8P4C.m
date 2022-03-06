clear;
% Diffusion coefficient
kappa = 0.02;

% Exact solution
%u_exact = @(x,t) sin(2*pi*x).*exp(-4*pi*pi*t*kappa);
u_exact = @(x,t) erfc(x/(sqrt(4*kappa*t)));

% Generate boundary condition and initial condition from exact solution
u_init = @(x,N) u_in(x,N);
gneg1 = @(t) u_exact(-1,t); 
g1 = @(t) u_exact(1,t);

% Spatial parameters
% Run over a sequence of grids to do convergence analysis
NN = [16; 32; 64; 128];

% Data structures for our error
L1_err = 0*NN;
L2_err = 0*NN;
max_err = 0*NN;

% Final time to run
T = 1.0;

% Run over all the grid sizes, compute error for each grid.
for i = 1:length(NN)
    
    % Generate grid
    x = linspace(-1,1,NN(i))';
    h = x(2) - x(1);
    
    % Timestep
    dt = 0.5 * h;
    
    % Run the diffusion solver
    U = solve_diff_tr(kappa, gneg1, g1, u_init, NN(i), T, dt);
    
    % Calculate error
    L1_err(i) = sum(abs(U(:,end) - u_exact(x,T))) * h;
    L2_err(i) = sqrt(sum(h*(U(:,end) - u_exact(x,T)).^2));
    max_err(i) = max(max(abs(U(:,end) - u_exact(x,T))));
end

% Plot the approximate and exact solution
figure(1); clf; hold on;
plot(x, U(:,1:4:end), 'k', 'linewidth', 3)
plot(x, u_exact(x,T), '--r', 'linewidth', 3)

% Make a convergence plot
figure(2); clf; hold on;
plot(1.0./NN, L1_err, 'linewidth', 4);
plot(1.0./NN, L2_err, 'linewidth', 4);
plot(1.0./NN, max_err, 'linewidth', 4);



% Calculate approximate orders by fitting a line to the errors.
l1_order = polyfit(log(1.0./NN), log(L1_err), 1);
l2_order = polyfit(log(1.0./NN), log(L2_err), 1);
max_order = polyfit(log(1.0./NN), log(max_err), 1);
l = legend("L1 order: " + num2str(l1_order(1)), "L2 order: " + num2str(l2_order(1)), "max order: " + num2str(max_order(1)));
set(l, 'fontsize', 18, 'location', 'southeast');
set(gca, 'xscale', 'log', 'yscale', 'log', 'fontsize', 18)
title('Convergence', 'fontsize', 18);


function U = solve_diff_tr(kappa, g0, g1, u_init, N, T, dt)
%%% Solve the diffusion equation using a trapezoidal rule

% Generate grid
x = linspace(-1,1,N)';
h = x(2) - x(1);
times = 0:dt:T;

% Initialize solution data structure.
% Each column is a solution at a particular time.
U = zeros(N, length(times));
U(:,1) = u_init(x,N);

% Note we can form the matrices outside the for loop since they don't
% change over time.
% Form A matrix
A = kappa/h^2*(-2*diag(ones(N-2,1),0) + diag(ones(N-3,1),1) + diag(ones(N-3,1),-1));
% Form matrix that is used to solve the system
mat = eye(N-2,N-2) - 0.5 * dt * A;
for i = 1:(length(times)-1)
    t = times(i);
    % g vector. This includes a fix up for the boundary conditions.
    gn = kappa/h^2*[g0(t); zeros(N-4,1); g1(t)];
    gn2 = kappa/h^2*[g0(t+dt); zeros(N-4,1); g1(t+dt)];
    % Form rhs for trapezoidal rule
    rhs = (eye(N-2,N-2) + 0.5 * dt * A)*U(2:end-1,i) + 0.5 * dt * (gn + gn2);
    % Trapezoidal Rule
    U(2:end-1,i+1) = mat\rhs;
    % Modify solution vector for the boundary conditions
    U(1,i+1) = g0(t+dt);
    U(end,i+1) = g1(t+dt);
end
end

function U = solve_diff_be(kappa, g0, g1, u_init, N, T, dt)
%%% Use backward Euler to solve the diffusion equation
x = linspace(-1,1,N)';
h = x(2) - x(1);
times = 0:dt:T;

U = zeros(N, length(times));
U(:,1) = u_init(x,N);

% Form A matrix
A = kappa/h^2*(-2*diag(ones(N-2,1),0) + diag(ones(N-3,1),1) + diag(ones(N-3,1),-1));
% Form solve matrix
mat = eye(N-2,N-2) - dt * A;
for i = 1:(length(times)-1)
    t = times(i);
    % Form g vector
    gn = kappa/h^2*[g0(t+dt); zeros(N-4,1); g1(t+dt)];
    % Form rhs
    rhs = (eye(N-2,N-2))*U(2:end-1,i) + dt * (gn);
    U(2:end-1,i+1) = linsolve(mat, rhs);%mat\rhs;
    U(1,i+1) = g0(t+dt);
    U(end,i+1) = g1(t+dt);
end
end


function initialU = u_in(x,N)

initialU = zeros(N,1);

for i = 1: N
    
if x(i,1) < 0
    
    initialU(i,1) = 2;
else
    
    initialU(i,1) = 0;
  
end

end

return

end
