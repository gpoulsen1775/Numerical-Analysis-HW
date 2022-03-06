%NOTE: "I've put line breaks in terms of percent signs to demonstrate where
%I've inserted code to modify the assingment." -Grant


%x Boundaries
xlow = 0;
xup = 1;

%Parameter a
a = 1;

%Interior
m = 64;

%Total With Boundaries
N = m + 2;

%Create x & dx
x = linspace(xlow, xup, N)';
dx = x(2) - x(1);

%Time Boundaries
t = 0;
T = 8;

%Step Vector & Drawing frequency
dt = 0.8 * dx;
n_steps = round(T / dt);
n_draw = round(0.01*n_steps);

%u_init = @(x) sin(pi*x).^4;
%u_init = @(x) heaviside(x - 0.5);
u_init = @(x) exp(-100*(x-0.5).^2);
%u_init = @(x) exp(-20 * (x-2.0).^2) + exp(-(x-5.0).^2);

% Note the degrees of freedom of u.
% u(2:N) are degrees of freedom
% u(1) and u(N+1) are "ghost cell" values filled from periodic boundaries.


%Assortment of Solutions
u_lf = [u_init(x); 0]; % lf -- Lax-Friedrichs
u_lw = [u_init(x); 0]; % lw -- Lax-Wendroff
u_uw = [u_init(x); 0]; % uw -- Upwind
u_lp = [u_init(x); 0]; % lp -- leapfrog
u_lp2 = u_lp;



%Modifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We Implement BW W/ a>0
u_bw = [u_init(x); 0; 0];

%Set Ghost Cells to Initial
u_bw(length(u_bw)-1) = u_bw(1);
u_bw(length(u_bw)) = u_bw(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Periodic boundaries
u_lf = fill_boundaries(u_lf);
u_lw = fill_boundaries(u_lw);
u_uw = fill_boundaries(u_uw);
u_lp = fill_boundaries(u_lp);


% Form index helper variable
I = 2:N;



%Modifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%We Make Helper to only alter
%non-ghost cells
I2 = 3:N+2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%



figure(1); clf;

for n = 1:n_steps
    
    % Perform update
    u_lf(I) = u_lf(I) - 0.5 * dt * a / dx * (u_lf(I + 1) - u_lf(I - 1)) + 0.5 * (u_lf(I+1) - 2*u_lf(I) + u_lf(I-1));
    u_lw(I) = u_lw(I) - 0.5 * dt * a / dx * (u_lw(I + 1) - u_lw(I - 1)) + 0.5 * a^2 * dt^2 / dx^2 * (u_lw(I+1) - 2*u_lw(I) + u_lw(I-1));
    u_uw(I) = u_uw(I) - dt * a / dx * (u_uw(I) - u_uw(I - 1));
        
    
    
    %Modifications
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Implement BW Iteration
    u_bw(I2) = u_bw(I2) - 0.5 * a * dt / dx * (3 * u_bw(I2) -4 * u_bw(I2-1) + u_bw(I2-2)) + (0.5 * a^2 * dt^2 / dx^2) * ( u_bw(I2) - 2 * u_bw(I2-1) + u_bw(I2-2));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    if (n == 1)
        
        % Note we need first timestep for leapfrog, use upwind
        u_lp2 = u_lp;
        u_lp(I) = u_lp(I) - dt * a / dx * (u_lp(I) - u_lp(I - 1));
        
    else
        % Update for leapfrog
        u_lp_new = u_lp2(I) - dt * a / dx * (u_lp(I+1) - u_lp(I - 1));
        u_lp2 = u_lp;
        u_lp(I) = u_lp_new;
        
    end
    
    % Fill boundary conditions
    u_lf = fill_boundaries(u_lf);
    u_lw = fill_boundaries(u_lw);
    u_uw = fill_boundaries(u_uw);
    u_lp = fill_boundaries(u_lp);
    
    
    
    %Modifications
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Reset Ghost Values 
    u_bw(1) = u_bw(length(u_bw)-1);
    u_bw(2) = u_bw(length(u_bw));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    % Update time
    t = t + dt;
    
    if (mod(n, n_draw) == 0 || n == n_steps)
        
        clf; hold on;
        
        % Draw solutions
        plot(x, u_lf(1:end-1), 'linewidth', 4);
        plot(x, u_lw(1:end-1), 'k', 'linewidth', 4);
        plot(x, u_uw(1:end-1), 'g', 'linewidth', 4);
        plot(x, u_lp(1:end-1), 'm', 'linewidth', 4);
        
        plot(x,u_bw(1:end-2), 'c', 'linewidth', 4);
        
        % Plot exact
        plot(x, u_exact(u_init, x, a, t, xup), 'r--', 'linewidth', 4);
        title(['t = ', num2str(t)]);
        legend('Lax-Friedrich', 'Lax-Wendroff', 'Upwind', 'Leapfrog','Beam-Warming');
        
        set(gca, 'fontsize', 18);
        pause(0.1);
    end
end

function u = u_exact(u_init, x, a, t, xmax)

% Exact solution. Assumes xlow = 0
xat = rem(x - a*t, xmax);
ineg = find(xat < 0);
xat(ineg) = xat(ineg) + xmax;
u = u_init(xat);

end

function u = fill_boundaries(u)

u(1) = u(end-1);
u(end) = u(2);

end