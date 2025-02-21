% RIS Single-User Grouping Method Simulation

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

% Test locations (as per paper Fig. 7)
angles = [70 90 130 145];  % degrees
distances = [170];  % cm (fixed distance for grouping comparison)

% Group sizes to test
group_sizes = [1 2 4 8];  % 1 means no grouping (iterative method)

% Convert to radians
angles_rad = deg2rad(angles);

% Initialize arrays for results
received_power = zeros(length(angles), length(group_sizes));
optimization_iterations = zeros(length(group_sizes), 1);

% Fixed Tx position (as per paper: 78°, 100cm)
tx_angle = deg2rad(78);
tx_distance = 100;
tx_pos = [tx_distance*cos(tx_angle), tx_distance*sin(tx_angle), 0];

% Run simulation for all test positions and group sizes
for a = 1:length(angles)
    % Calculate Rx position
    rx_pos = [distances*cos(angles_rad(a)), ...
             distances*sin(angles_rad(a)), 0];
    
    % Calculate channels
    h = get_tx_ris_channel(N, Nx, Ny, tx_pos, dx, dy, k);
    g = get_ris_rx_channel(N, Nx, Ny, rx_pos, dx, dy, k);
    
    % Test different group sizes
    for g_idx = 1:length(group_sizes)
        G = group_sizes(g_idx);
        [theta_opt, max_power, num_iter] = group_optimization(h, g, G, Nx, Ny);
        received_power(a, g_idx) = 10*log10(max_power);
        optimization_iterations(g_idx) = num_iter;
    end
end

% Plot results
figure;
for a = 1:length(angles)
    plot(1:length(group_sizes), received_power(a,:), '-o', 'LineWidth', 1.5);
    hold on;
end
hold off;
xlabel('Group Size Index');
xticks(1:length(group_sizes));
xticklabels(arrayfun(@num2str, group_sizes, 'UniformOutput', false));
ylabel('Received Power (dB)');
title('RIS Performance with Different Grouping Sizes');
legend(arrayfun(@(x) [num2str(x) '°'], angles, 'UniformOutput', false));
grid on;

% Plot number of iterations
figure;
bar(optimization_iterations);
xlabel('Group Size');
xticklabels(arrayfun(@num2str, group_sizes, 'UniformOutput', false));
ylabel('Number of Iterations');
title('Number of Optimization Iterations vs Group Size');
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

function [theta_opt, max_power, num_iter] = group_optimization(h, g, G, Nx, Ny)
    N = length(h);
    theta_opt = zeros(N, 1);  % Initialize phases
    max_power = 0;
    
    % Calculate number of groups
    num_groups = ceil(N/G);
    num_iter = num_groups * 2;  % 2 phase states per group
    
    % Create group indices
    group_indices = cell(num_groups, 1);
    current_idx = 1;
    
    % Assign elements to groups based on spatial proximity
    for group = 1:num_groups
        if group == num_groups
            group_indices{group} = current_idx:N;  % Last group gets remaining elements
        else
            group_indices{group} = current_idx:min(current_idx+G-1, N);
        end
        current_idx = current_idx + G;
    end
    
    % Optimize phase for each group
    for group = 1:num_groups
        max_local_power = -inf;
        best_phase = 0;
        
        % Try both phase states (0° and 180°)
        for phase = [0 pi]
            % Apply same phase to all elements in the group
            theta_opt(group_indices{group}) = phase;
            Theta = diag(exp(1j*theta_opt));
            received_power = abs(g'*Theta*h)^2;
            
            if received_power > max_local_power
                max_local_power = received_power;
                best_phase = phase;
            end
        end
        
        % Apply best phase to all elements in the group
        theta_opt(group_indices{group}) = best_phase;
        max_power = max_local_power;
    end
end