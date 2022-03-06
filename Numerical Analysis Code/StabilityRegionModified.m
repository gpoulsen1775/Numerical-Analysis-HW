clear;

% We've Modified the Code to Construct the region of stability for the
% Linear Multistep Method Given in Problem 4. We Consider Different Values
% of Theta...

% Theta of Choice
theta = 1; 

% Stability Argument
R = @(z) (1 + z -z*theta) ./ (1 - z*theta); 

% Bounds
bounds = [-4, 2, -3, 3];

% Function For Plotting Region of Stability
plotStabilityRegion(R, bounds);

% NOTE: THE CODE SWITCHES WHAT PORTION OF THE MAP IS COLORED, SO THE FIGURE
% IS MORESO AMBIGIOUS BUT THE REGION OF STABILITY IS ALWAYS THE SECTION
% SURROUNDED BY THE NEON GREEN BOUNDARY!!!

function plotStabilityRegion(R, bounds)

% This only works for one step methods!
xlow = bounds(1);
xup = bounds(2);
ylow = bounds(3);
yup = bounds(4);

[X, Y] = meshgrid(linspace(xlow, xup, 500), linspace(ylow, yup, 500));
Z = X + 1i * Y;
Rvals = abs(R(Z));

figure(1); clf; hold on;

surf(X,Y, Rvals, 'Edgecolor', 'none');
colormap([0 0 0; 1 1 1]) % Create a two color colormap.
caxis([.99 1.01]) % Adjust colormap so we only plot colors between .99 and 1.01
contour(X,Y,Rvals,[1 1],'g', 'linewidth', 4) % Plot the contour at 1
view(0, -90); % View from the bottom so we can see the contour
title('Stability Region', 'fontsize', 18);
xlabel('Real')
ylabel('Imaginary')
set(gca, 'fontsize', 18);

end