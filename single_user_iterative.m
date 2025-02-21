% RIS Single-User Iterative Method Simulation

% System Parameters
fc = 5.2e9;                  % Carrier frequency (5.2 GHz)
c = 3e8;                     % Speed of light
lambda = c/fc;               % Wavelength
k = 2*pi/lambda;            % Wave number

% RIS Parameters
Nx = 10;                     % Number of elements in x-direction
Ny = 8;                      % Number of elements in y-direction
N = 76;                      % Total number of elements (removing 2x2 corner)
dx = lambda/2;              % Element spacing in x-direction
dy = lambda/2;              % Element spacing in y-direction

% Test locations
angles = [50 70 90 110 130 145];  % degrees
distances = [120 170 220 270 320 420];  % cm

% Convert to radians
angles_rad = deg2rad(angles);

% Initialize arrays for results
received_power = zeros(length(angles), length(distances));
optimized_power = zeros(length(angles), length(distances));

% Fixed Tx position (as per paper: 78°, 100cm)
tx_angle = deg2rad(78);
tx_distance = 100;
tx_pos = [tx_distance*cos(tx_angle), tx_distance*sin(tx_angle), 0];

% Run simulation for all test positions
for a = 1:length(angles)
    for d = 1:length(distances)
        % Calculate Rx position
        rx_pos = [distances(d)*cos(angles_rad(a)), ...
                 distances(d)*sin(angles_rad(a)), 0];
        
        % Calculate channels
        h = get_tx_ris_channel(N, Nx, Ny, tx_pos, dx, dy, k);
        g = get_ris_rx_channel(N, Nx, Ny, rx_pos, dx, dy, k);
        
        % Calculate received power without RIS optimization
        Theta_init = eye(N);
        received_power(a,d) = 10*log10(abs(g'*Theta_init*h)^2);
        
        % Optimize RIS phases
        [theta_opt, max_power] = iterative_optimization(h, g);
        optimized_power(a,d) = 10*log10(max_power);
    end
end

% Plot results
figure;
for a = 1:length(angles)
    subplot(2,3,a);
    plot(distances, received_power(a,:), 'b--', 'LineWidth', 1.5);
    hold on;
    plot(distances, optimized_power(a,:), 'r-', 'LineWidth', 1.5);
    title(['Angle = ' num2str(angles(a)) '°']);
    xlabel('Distance (cm)');
    ylabel('Received Power (dB)');
    legend('Without Optimization', 'With Optimization');
    grid on;
end

sgtitle('RIS Performance with Iterative Optimization');

% Local function definitions below
function h = get_tx_ris_channel(N, Nx, Ny, tx_pos, dx, dy, k)
    h = zeros(N, 1);
    idx = 1;
    for nx = 1:Nx
        for ny = 1:Ny
            if ~(nx >= 9 && ny >= 7)  % Skip the 2x2 corner
                % Calculate element position
                x = (nx-1)*dx;
                y = (ny-1)*dy;
                
                % Calculate distance and phase
                d = sqrt((tx_pos(1)-x)^2 + (tx_pos(2)-y)^2 + tx_pos(3)^2);
                phase = exp(-1j*k*d);
                
                % Path loss
                pl = sqrt(1/(4*pi*d^2));
                
                h(idx) = pl * phase;
                idx = idx + 1;
            end
        end
    end
end

function g = get_ris_rx_channel(N, Nx, Ny, rx_pos, dx, dy, k)
    g = zeros(N, 1);
    idx = 1;
    for nx = 1:Nx
        for ny = 1:Ny
            if ~(nx >= 9 && ny >= 7)  % Skip the 2x2 corner
                % Calculate element position
                x = (nx-1)*dx;
                y = (ny-1)*dy;
                
                % Calculate distance and phase
                d = sqrt((rx_pos(1)-x)^2 + (rx_pos(2)-y)^2 + rx_pos(3)^2);
                phase = exp(-1j*k*d);
                
                % Path loss
                pl = sqrt(1/(4*pi*d^2));
                
                g(idx) = pl * phase;
                idx = idx + 1;
            end
        end
    end
end

function [theta_opt, max_power] = iterative_optimization(h, g)
    N = length(h);
    theta_opt = zeros(N, 1);  % Initialize phases
    max_power = 0;
    
    for n = 1:N
        max_local_power = -inf;
        best_phase = 0;
        
        % Try both phase states (0° and 180°)
        for phase = [0 pi]
            theta_opt(n) = phase;
            Theta = diag(exp(1j*theta_opt));
            received_power = abs(g'*Theta*h)^2;
            
            if received_power > max_local_power
                max_local_power = received_power;
                best_phase = phase;
            end
        end
        
        theta_opt(n) = best_phase;
        max_power = max_local_power;
    end
end