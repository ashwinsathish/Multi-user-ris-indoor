% RIS Multi-User Physical Layer Security Simulation

% System Parameters
fc = 5.2e9;                  % Carrier frequency (5.2 GHz)
c = 3e8;                     % Speed of light
lambda = c/fc;               % Wavelength
k = 2*pi/lambda;            % Wave number

% RIS Parameters
Nx = 10;                     % Number of elements in x-direction
Ny = 8;                      % Number of elements in y-direction
N = 76;                      % Total elements (removing 2x2 corner)
N_half = floor(N/2);        % Elements per user after partitioning
dx = lambda/2;              % Element spacing in x-direction
dy = lambda/2;              % Element spacing in y-direction

fprintf('Total elements: %d, Elements per half: %d\n', N, N_half);

% Codebook reference locations (at 170cm as per paper)
codebook_angles = [50 70 90 110 130 145];  % degrees
codebook_distance = 170;  % cm

% Movement Scenario 1 (Fig. 12 - similar patterns)
rx1_angles_s1 = [90 90 90 90 90 90 90];    % Rx1 angles
rx1_dist_s1 = [120 170 220 270 320 370 420]; % Rx1 distances
rx2_angles_s1 = [110 110 110 110 110 110 110]; % Rx2 angles
rx2_dist_s1 = [120 170 220 270 320 370 420];   % Rx2 distances

% Movement Scenario 2 (Fig. 13 - opposite directions)
rx1_angles_s2 = [50 70 90 110 130 145 130 110 90 70];  % Rx1 angles
rx1_dist_s2 = [170 170 170 170 170 170 220 220 220 220];  % Rx1 distances
rx2_angles_s2 = [145 130 110 90 70 50 70 90 110 130];    % Rx2 angles
rx2_dist_s2 = [220 220 220 220 220 220 170 170 170 170];  % Rx2 distances

% Fixed Tx position (as per paper: 78°, 105cm)
tx_angle = deg2rad(78);
tx_distance = 105;
tx_pos = [tx_distance*cos(tx_angle), tx_distance*sin(tx_angle), 0];

% Generate codebooks for both halves of RIS
fprintf('Generating codebooks...\n');
codebook_rx1 = cell(length(codebook_angles), 5);  % 5 configurations per location
codebook_rx2 = cell(length(codebook_angles), 5);
codebook_powers_rx1 = zeros(length(codebook_angles), 5);
codebook_powers_rx2 = zeros(length(codebook_angles), 5);

for a = 1:length(codebook_angles)
    % Calculate positions
    rx_pos = [codebook_distance*cos(deg2rad(codebook_angles(a))), ...
              codebook_distance*sin(deg2rad(codebook_angles(a))), 0];
    
    % Calculate channels for both halves
    [h1, h2] = get_split_tx_ris_channel(N_half, Nx, Ny, tx_pos, dx, dy, k);
    [g1, g2] = get_split_ris_rx_channel(N_half, Nx, Ny, rx_pos, dx, dy, k);
    
    % Generate 5 different configurations
    for i = 1:5
        % Randomize initial phases to get different configurations
        [theta1, power1] = iterative_optimization(h1, g1, rand(N_half,1)*pi);
        [theta2, power2] = iterative_optimization(h2, g2, rand(N_half,1)*pi);
        
        codebook_rx1{a,i} = theta1;
        codebook_rx2{a,i} = theta2;
        codebook_powers_rx1(a,i) = power1;
        codebook_powers_rx2(a,i) = power2;
    end
    fprintf('Codebook entry %d/%d completed\n', a, length(codebook_angles));
end

% Run simulations for both scenarios
[sc1, rx1_p1, rx2_p1] = simulate_scenario(rx1_angles_s1, rx1_dist_s1, ...
    rx2_angles_s1, rx2_dist_s1, N_half, Nx, Ny, tx_pos, dx, dy, k, ...
    codebook_angles, codebook_rx1, codebook_rx2);

[sc2, rx1_p2, rx2_p2] = simulate_scenario(rx1_angles_s2, rx1_dist_s2, ...
    rx2_angles_s2, rx2_dist_s2, N_half, Nx, Ny, tx_pos, dx, dy, k, ...
    codebook_angles, codebook_rx1, codebook_rx2);

% Plot results for Scenario 1
figure;
subplot(3,1,1);
plot(1:length(rx1_angles_s1), rx1_p1.ris_off, 'b--', 'DisplayName', 'RIS-off');
hold on;
plot(1:length(rx1_angles_s1), rx1_p1.codebook, 'r-', 'DisplayName', 'Codebook');
plot(1:length(rx1_angles_s1), rx1_p1.iterative, 'g-', 'DisplayName', 'Iterative');
title('Scenario 1: Rx1 (Intended User) Power');
ylabel('Received Power (dB)');
legend('Location', 'best');
grid on;

subplot(3,1,2);
plot(1:length(rx1_angles_s1), rx2_p1.ris_off, 'b--', 'DisplayName', 'RIS-off');
hold on;
plot(1:length(rx1_angles_s1), rx2_p1.codebook, 'r-', 'DisplayName', 'Codebook');
plot(1:length(rx1_angles_s1), rx2_p1.iterative, 'g-', 'DisplayName', 'Iterative');
title('Rx2 (Unintended User) Power');
ylabel('Received Power (dB)');
legend('Location', 'best');
grid on;

subplot(3,1,3);
plot(1:length(rx1_angles_s1), sc1.ris_off, 'b--', 'DisplayName', 'RIS-off');
hold on;
plot(1:length(rx1_angles_s1), sc1.codebook, 'r-', 'DisplayName', 'Codebook');
plot(1:length(rx1_angles_s1), sc1.iterative, 'g-', 'DisplayName', 'Iterative');
title('Secrecy Capacity');
xlabel('Position Index');
ylabel('Secrecy Capacity (bps/Hz)');
legend('Location', 'best');
grid on;

% Plot results for Scenario 2
figure;
subplot(3,1,1);
plot(1:length(rx1_angles_s2), rx1_p2.ris_off, 'b--', 'DisplayName', 'RIS-off');
hold on;
plot(1:length(rx1_angles_s2), rx1_p2.codebook, 'r-', 'DisplayName', 'Codebook');
plot(1:length(rx1_angles_s2), rx1_p2.iterative, 'g-', 'DisplayName', 'Iterative');
title('Scenario 2: Rx1 (Intended User) Power');
ylabel('Received Power (dB)');
legend('Location', 'best');
grid on;

subplot(3,1,2);
plot(1:length(rx1_angles_s2), rx2_p2.ris_off, 'b--', 'DisplayName', 'RIS-off');
hold on;
plot(1:length(rx1_angles_s2), rx2_p2.codebook, 'r-', 'DisplayName', 'Codebook');
plot(1:length(rx1_angles_s2), rx2_p2.iterative, 'g-', 'DisplayName', 'Iterative');
title('Rx2 (Unintended User) Power');
ylabel('Received Power (dB)');
legend('Location', 'best');
grid on;

subplot(3,1,3);
plot(1:length(rx1_angles_s2), sc2.ris_off, 'b--', 'DisplayName', 'RIS-off');
hold on;
plot(1:length(rx1_angles_s2), sc2.codebook, 'r-', 'DisplayName', 'Codebook');
plot(1:length(rx1_angles_s2), sc2.iterative, 'g-', 'DisplayName', 'Iterative');
title('Secrecy Capacity');
xlabel('Position Index');
ylabel('Secrecy Capacity (bps/Hz)');
legend('Location', 'best');
grid on;

% Function definitions
function [secrecy_capacity, rx1_power, rx2_power] = simulate_scenario(rx1_angles, rx1_dist, ...
    rx2_angles, rx2_dist, N_half, Nx, Ny, tx_pos, dx, dy, k, ...
    codebook_angles, codebook_rx1, codebook_rx2)
    
    num_points = length(rx1_angles);
    secrecy_capacity = struct('ris_off', zeros(num_points,1), ...
                            'codebook', zeros(num_points,1), ...
                            'iterative', zeros(num_points,1));
    rx1_power = struct('ris_off', zeros(num_points,1), ...
                     'codebook', zeros(num_points,1), ...
                     'iterative', zeros(num_points,1));
    rx2_power = struct('ris_off', zeros(num_points,1), ...
                     'codebook', zeros(num_points,1), ...
                     'iterative', zeros(num_points,1));
                     
    for p = 1:num_points
        % Calculate positions
        rx1_pos = [rx1_dist(p)*cos(deg2rad(rx1_angles(p))), ...
                  rx1_dist(p)*sin(deg2rad(rx1_angles(p))), 0];
        rx2_pos = [rx2_dist(p)*cos(deg2rad(rx2_angles(p))), ...
                  rx2_dist(p)*sin(deg2rad(rx2_angles(p))), 0];
        
        % Calculate channels
        [h1, h2] = get_split_tx_ris_channel(N_half, Nx, Ny, tx_pos, dx, dy, k);
        [g1_rx1, g2_rx1] = get_split_ris_rx_channel(N_half, Nx, Ny, rx1_pos, dx, dy, k);
        [g1_rx2, g2_rx2] = get_split_ris_rx_channel(N_half, Nx, Ny, rx2_pos, dx, dy, k);
        
        % 1. RIS-off performance
        theta_off = zeros(N_half, 1);
        rx1_power.ris_off(p) = calculate_received_power(h1, g1_rx1, theta_off) + ...
                               calculate_received_power(h2, g2_rx1, theta_off);
        rx2_power.ris_off(p) = calculate_received_power(h1, g1_rx2, theta_off) + ...
                               calculate_received_power(h2, g2_rx2, theta_off);
        secrecy_capacity.ris_off(p) = calculate_secrecy_capacity(rx1_power.ris_off(p), ...
                                                               rx2_power.ris_off(p));
        
        % 2. Codebook performance
        [~, idx1] = min(abs(codebook_angles - rx1_angles(p)));
        [~, idx2] = min(abs(codebook_angles - rx2_angles(p)));
        
        max_secrecy = -inf;
        best_power_rx1 = 0;
        best_power_rx2 = 0;
        
        for i = 1:5
            for j = 1:5
                theta1 = codebook_rx1{idx1,i};
                theta2 = codebook_rx2{idx2,j};
                
                p_rx1 = calculate_received_power(h1, g1_rx1, theta1) + ...
                        calculate_received_power(h2, g2_rx1, theta2);
                p_rx2 = calculate_received_power(h1, g1_rx2, theta1) + ...
                        calculate_received_power(h2, g2_rx2, theta2);
                
                sec_cap = calculate_secrecy_capacity(p_rx1, p_rx2);
                
                if sec_cap > max_secrecy
                    max_secrecy = sec_cap;
                    best_power_rx1 = p_rx1;
                    best_power_rx2 = p_rx2;
                end
            end
        end
        
        rx1_power.codebook(p) = best_power_rx1;
        rx2_power.codebook(p) = best_power_rx2;
        secrecy_capacity.codebook(p) = max_secrecy;
        
        % 3. Online iterative optimization
        [theta1, theta2, p_rx1, p_rx2] = joint_iterative_optimization(h1, h2, ...
            g1_rx1, g2_rx1, g1_rx2, g2_rx2);
        
        rx1_power.iterative(p) = p_rx1;
        rx2_power.iterative(p) = p_rx2;
        secrecy_capacity.iterative(p) = calculate_secrecy_capacity(p_rx1, p_rx2);
    end
end

function [h1, h2] = get_split_tx_ris_channel(N_half, Nx, Ny, tx_pos, dx, dy, k)
    h1 = zeros(N_half, 1);
    h2 = zeros(N_half, 1);
    idx1 = 1;
    idx2 = 1;
    
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
                
                % Assign to first or second half
                if idx1 <= N_half
                    h1(idx1) = pl * phase;
                    idx1 = idx1 + 1;
                elseif idx2 < N_half
                    h2(idx2) = pl * phase;
                    idx2 = idx2 + 1;
                end
            end
        end
    end
end

function [g1, g2] = get_split_ris_rx_channel(N_half, Nx, Ny, rx_pos, dx, dy, k)
    g1 = zeros(N_half, 1);
    g2 = zeros(N_half, 1);
    idx1 = 1;
    idx2 = 1;
    
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
                
                % Assign to first or second half
                if idx1 <= N_half
                    g1(idx1) = pl * phase;
                    idx1 = idx1 + 1;
                elseif idx2 < N_half
                    g2(idx2) = pl * phase;
                    idx2 = idx2 + 1;
                end
            end
        end
    end
end

function power = calculate_received_power(h, g, theta)
    % Ensure proper dimensions
    assert(length(h) == length(theta), 'h and theta must have same length');
    assert(length(g) == length(theta), 'g and theta must have same length');
    
    % Create diagonal matrix of phase shifts
    Theta = diag(exp(1j*theta));
    
    % Calculate received power
    power = 10*log10(abs(g'*Theta*h)^2);
end

function secrecy_capacity = calculate_secrecy_capacity(p_rx1, p_rx2)
    % Using equation (6) from the paper
    secrecy_capacity = 0.1*log2(10)*(p_rx1 - p_rx2);
end

function [theta_opt, max_power] = iterative_optimization(h, g, init_phase)
    N = length(h);
    theta_opt = init_phase;  % Use provided initial phases
    max_power = -inf;
    
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
        if max_local_power > max_power
            max_power = max_local_power;
        end
    end
    max_power = 10*log10(max_power);  % Convert to dB
end

function [theta1_opt, theta2_opt, p_rx1, p_rx2] = joint_iterative_optimization(...
    h1, h2, g1_rx1, g2_rx1, g1_rx2, g2_rx2)
    
    % Initialize phases randomly
    theta1_opt = rand(length(h1), 1) * pi;
    theta2_opt = rand(length(h2), 1) * pi;
    
    % Iterative optimization for maximum secrecy capacity
    max_secrecy = -inf;
    max_iterations = 10;  % Limit iterations for practical implementation
    p_rx1 = 0;
    p_rx2 = 0;
    
    for iter = 1:max_iterations
        % Optimize first half (for Rx1)
        [theta1_new, ~] = iterative_optimization(h1, g1_rx1, theta1_opt);
        
        % Optimize second half (against Rx2)
        [theta2_new, ~] = iterative_optimization(h2, g2_rx2, theta2_opt);
        
        % Calculate total received powers
        p_rx1_new = calculate_received_power(h1, g1_rx1, theta1_new) + ...
                    calculate_received_power(h2, g2_rx1, theta2_new);
        p_rx2_new = calculate_received_power(h1, g1_rx2, theta1_new) + ...
                    calculate_received_power(h2, g2_rx2, theta2_new);
        
        secrecy = calculate_secrecy_capacity(p_rx1_new, p_rx2_new);
        
        if secrecy > max_secrecy
            max_secrecy = secrecy;
            theta1_opt = theta1_new;
            theta2_opt = theta2_new;
            p_rx1 = p_rx1_new;
            p_rx2 = p_rx2_new;
        end
    end
end