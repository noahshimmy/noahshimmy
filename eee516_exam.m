%% Runge-Kutta method 
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

%% GAUSSIAN ELIMINATION 
clc % clear the command window for a fresh output

% --- Defining the System (Ax = b) ---
A = [6 -1; -2 7]; % Coefficient matrix A (represents the left side of equations)
b = [10; 50];     % Constant vector b (represents the right side of equations)
n = length(b);   % Find the number of equations/unknown variables

% --- Step 1: Forward Elimination ---
% The goal here is to make 'A' an upper triangular matrix (zeros below the main diagonal)
for k = 1:n-1 % Loop through the pivot equations (columns)
    for i = k+1:n % Loop through the equations below the pivot (rows)
        
        % Calculate the multiplier needed to eliminate the variable in the current row
        factor = A(i,k)/A(k,k); 
        
        % Perform the row operation: R_i = R_i - factor * R_k
        A(i,k:n) = A(i,k:n) - factor*A(k,k:n); % Update matrix A
        b(i) = b(i) - factor*b(k);             % Update vector b to match
        
    end
end

% --- Step 2: Back Substitution ---
% Now that we have a triangular matrix, we solve from the bottom equation up to the top
v = zeros(n,1);         % Pre-allocate the solution vector 'v' with zeros
v(n) = b(n)/A(n,n);     % Solve the very last equation directly (since it only has one unknown left)

for i = n-1:-1:1 % Loop backwards from the second-to-last equation up to the first equation
    
    % Substitute the already-found 'v' values into the current row to solve for the unknown 'v(i)'
    v(i) = (b(i) - A(i,i+1:n)*v(i+1:n))/A(i,i); 
    
end

% --- Display results ---
disp('Solutions:'); % Print text header
disp(v);            % Print the final calculated values for the variables

%% EULER METHOD
clc       % clear the command window
clear all % clear all stored variables from the workspace

% --- Initial Conditions and Parameters ---
y0 = 1;   % the initial value of y when t = 0
h = 0.2;  % step size (time increment)
t = 0:h:2; % create a time vector from 0 to 2, stepping by h

% --- Exact Solution (for comparison) ---
% The analytical solution to the ODE dy/dt = t - 2y
yexact = 0.25*(2*t - 1 + 5*exp(-2*t)); 

% --- Pre-allocating Vectors ---
ystar = zeros(size(t)); % create an array of zeros the same size as 't' to store our answers
ystar(1) = y0;          % FIX: set ONLY the first element to the initial condition

% --- Euler Method Main Loop ---
for i = 1:(length(t)-1) % loop through each time step
    
    % Calculate the slope at the current point: f(t, y) = t - 2y
    k1 = t(i) - 2*ystar(i); 
    
    % Predict the next y value using the Euler formula: y_new = y_old + slope * step_size
    ystar(i+1) = ystar(i) + k1*h; 
    
end

% --- Plotting the Results ---
plot(t, yexact, t, ystar, '--'); % plot exact as a solid line, approximate as dashed ('--')
legend('Exact', 'Approximate');  % add a legend so we know which line is which

% --- Displaying the Data Table ---
disp('     time   y_approx   yexact'); % print column headers
disp([t; ystar; yexact]');           % combine the three row vectors, transpose ('), and print as columns

 %% Y-Bus Formation
clc
% Line impedances
Z12 = 0.02 + 0.06i;
Z23 = 0.08 + 0.24i;
Z13 = 0.06 + 0.18i;
% Admittances
Y12 = 1/Z12; Y23 = 1/Z23; Y13 = 1/Z13;
% Y-bus formation
Y = zeros(3,3);
Y(1,1) = Y12 + Y13;
Y(2,2) = Y12 + Y23;
Y(3,3) = Y13 + Y23;
Y(1,2) = -Y12; Y(2,1) = -Y12;
Y(2,3) = -Y23; Y(3,2) = -Y23;
Y(1,3) = -Y13; Y(3,1) = -Y13;
disp('Y-bus matrix:');
disp(Y)