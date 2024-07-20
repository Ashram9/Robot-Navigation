% All-epoch Positioning
% Define constants
Define_Constants

% Load .csv data file
filename = 'Workshop1_Pseudo_ranges.csv';
table_pseudo_ranges = readtable(filename);
% 从CSV文件中读取卫星编号
satellite_numbers = table_pseudo_ranges{1, 2:end};

% Get the number of epochs
num_epochs = size(table_pseudo_ranges, 1) - 1; % Subtract 1 for the header

% Initialize arrays to store results for all epochs
all_positions = zeros(num_epochs, 3);
all_clock_offsets = zeros(num_epochs, 1);

% Initialize initial position at Earth center
r_ea_e = [0; 0; 0];
% Initialize receiver clock offset
receiver_clock_offset = 0;

% Initialize tolerance for convergence
tolerance = 0.1; % 10cm tolerance

for epoch = 1:num_epochs
    % Extract pseudo-ranges for the current epoch
    time = epoch - 1; % Time starts from 0
    pseudo_ranges = table_pseudo_ranges{1 + epoch, 2:end};
    
    iteration = 1;
    while true
        % Predicted state vector
        x_prev = [r_ea_e; receiver_clock_offset];

        % Calculate satellite positions and other parameters
        r_es_e = zeros(3, numel(pseudo_ranges));
        r_as = zeros(1, numel(pseudo_ranges));
        u_as_e = zeros(3, numel(pseudo_ranges));
        delta_z = zeros(numel(pseudo_ranges), 1);
        H = zeros(numel(pseudo_ranges), 4);

        for j = 1:numel(pseudo_ranges)
            satellite_number = satellite_numbers(j);
            [r_ej_e, ~] = Satellite_position_and_velocity(time*60, satellite_number);
            r_es_e(:, j) = r_ej_e;
            
            C_I_e = eye(3);
            r_aj = norm(C_I_e * r_es_e(:, j) - r_ea_e);
            C_I_e = eye(3) - Omega_ie * r_aj / c;
            r_aj = norm(C_I_e * r_es_e(:, j) - r_ea_e);
            r_as(:, j) = r_aj;

            u_aj_e = (C_I_e * r_es_e(:, j) - r_ea_e) / r_aj;
            u_as_e(:, j) = u_aj_e;

            % Calculate measurements and measurement innovation
            delta_z(j) = pseudo_ranges(j) - (r_aj + receiver_clock_offset);
            u_aj_e = u_as_e(:, j);
            H(j, :) = [-u_aj_e', 1];
        end

        % Calculate Moore-Penrose pseudo-inverse of H
        H_pseudo_inv = pinv(H);

        % Calculate the new state vector using least-squares
        x_new = x_prev + H_pseudo_inv * delta_z;
        % Extract position and receiver clock offset
        r_ea_e = x_new(1:3);
        receiver_clock_offset = x_new(4);

        % Calculate the position difference between iterations
        position_difference = norm(r_ea_e - x_prev(1:3));
        % Check for convergence
        if position_difference < tolerance
            break;
        end
        iteration = iteration + 1;
    end

    % Store the results for this epoch
    all_positions(epoch, :) = r_ea_e';
    all_clock_offsets(epoch) = receiver_clock_offset;
end


% Display the results (you can also save them to a file)
disp('Time(s)   Latitude (deg)   Longitude (deg)   Height (m)');
for epoch = 1:num_epochs
    time = epoch - 1;
    position = all_positions(epoch, :);
    % Convert ECEF positions to latitude, longitude, and height
    [L_b, lambda_b, h_b, ~] = pv_ECEF_to_NED(position', [0; 0; 0]);
    % Convert latitude and longitude from radians to degrees
    L_b_deg = rad_to_deg * L_b;
    lambda_b_deg = rad_to_deg * lambda_b;
    disp([num2str(time*60), '       ', num2str(L_b_deg), '        ', ...
        num2str(lambda_b_deg), '       ', num2str(h_b)]);
end
