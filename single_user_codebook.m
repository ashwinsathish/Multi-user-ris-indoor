% RIS Single-User Codebook Method Simulation

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

% Codebook reference locations (at 170cm as per paper)
codebook_angles = [50 70 90 110 130 145];  % degrees
codebook_distance = 170;  % cm

% Test locations (interior points between reference locations)
test_angles = [60 80 100 120 135];  % degrees
test_distance = 170;  % cm

% Convert to radians
codebook_angles_rad = deg2rad(codebook_angles);
test_angles_rad = deg2rad(test_angles);

% Fixed Tx position (as per paper: 78°, 100cm)
tx_angle = deg2rad(78);
tx_distance = 100;
tx_pos = [tx_distance*cos(tx_angle), tx_distance*sin(tx_angle), 0];

% Initialize codebook
codebook = cell(length(codebook_angles), 1);
codebook_powers = zeros(length(codebook_angles), 1);

% Generate codebook
fprintf('Generating codebook...\n');
for a = 1:length(codebook_angles)
    % Calculate Rx position for codebook entry
    rx_pos = [codebook_distance*cos(codebook_angles_rad(a)), ...
             codebook_distance*sin(codebook_angles_rad(a)), 0];
    
    % Calculate channels
    h = get_tx_ris_channel(N, Nx, Ny, tx_pos, dx, dy, k);
    g = get_ris_rx_channel(N, Nx, Ny, rx_pos, dx, dy, k);
    
    % Generate optimal configuration using iterative method
    [theta_opt, max_power] = iterative_optimization(h, g);
    
    % Store in codebook
    codebook{a} = theta_opt;
    codebook_powers(a) = 10*log10(max_power);
    
    fprintf('Codebook entry %d/%d completed\n', a, length(codebook_angles));
end

% Initialize results arrays
ris_off_power = zeros(length(test_angles), 1);
codebook_power = zeros(length(test_angles), 1);
online_power = zeros(length(test_angles), 1);

% Test at interior points
fprintf('Testing at interior points...\n');
for a = 1:length(test_angles)
    % Calculate Rx test position
    rx_pos = [test_distance*cos(test_angles_rad(a)), ...
             test_distance*sin(test_angles_rad(a)), 0];
    
    % Calculate channels
    h = get_tx_ris_channel(N, Nx, Ny, tx_pos, dx, dy, k);
    g = get_ris_rx_channel(N, Nx, Ny, rx_pos, dx, dy, k);
    
    % 1. RIS-off performance
    Theta_off = eye(N);
    ris_off_power(a) = 10*log10(abs(g'*Theta_off*h)^2);
    
    % 2. Codebook performance - find nearest codebook entry
    [~, nearest_idx] = min(abs(codebook_angles - test_angles(a)));
    theta_codebook = codebook{nearest_idx};
    Theta_codebook = diag(exp(1j*theta_codebook));
    codebook_power(a) = 10*log10(abs(g'*Theta_codebook*h)^2);
    
    % 3. Online iterative performance
    [theta_online, max_power] = iterative_optimization(h, g);
    online_power(a) = 10*log10(max_power);
    
    fprintf('Test point %d/%d completed\n', a, length(test_angles));
end

% Plot results
figure;
plot(test_angles, ris_off_power, 'b--o', 'LineWidth', 1.5, 'DisplayName', 'RIS-off');
hold on;
plot(test_angles, codebook_power, 'r-s', 'LineWidth', 1.5, 'DisplayName', 'Codebook');
plot(test_angles, online_power, 'g-^', 'LineWidth', 1.5, 'DisplayName', 'Online Iterative');
plot(codebook_angles, codebook_powers, 'k*', 'LineWidth', 1.5, 'DisplayName', 'Codebook Reference');
hold off;

xlabel('Angle (degrees)');
ylabel('Received Power (dB)');
title('RIS Performance Comparison at 170cm');
legend('Location', 'best');
grid on;

% Local function definitions
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