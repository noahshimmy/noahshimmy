% Runge-Kutta method 
clc
v0 = 15;                % initial voltage/value
h = 0.1;                % step size
t = 0:h:1;              % time vector from 0 to 1 with step h

vexact = 15*exp(-5*t);  % analytical/exact solution for comparison
vstar = zeros(size(t)); % pre-allocate the vector with zeros
vstar(1) = v0;          % set the first element to the initial value (v0)

for i = 1:(length(t)-1) % loop through the time steps
    
    % Calculating the four RK4 increments
    k1 = h*(-5*vstar(i));
    k2 = h*(-5*(vstar(i) + 0.5*k1));
    k3 = h*(-5*(vstar(i) + 0.5*k2));
    k4 = h*(-5*(vstar(i) + k3));
    
    % Updating the next value in the array
    vstar(i+1) = vstar(i) + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
    
end % end of the loop

% Displaying final results
fprintf('Final RK4 Value: %.6f\n', vstar(end));
fprintf('Final Exact Value: %.6f\n', vexact(end));