function [r_ea_e, v_ea_e, receiver_clock_offset, receiver_clock_drift] = GNSS_LSE(k)

% GNSS LSE All-Epochs
addpath('Coursework 1 Data/');
addpath('Coursework 1 Software for Students to Use/');
% Define constants
Define_Constants

sigma_rho = 5; % Measurement error standard deviation
threshold = 6; % Outlier detection threshold

% Load .csv data file
table_pseudo_ranges = readtable('Pseudo_ranges.csv');
% Load .csv data files for pseudo-range rates
table_pseudo_range_rates = readtable('Pseudo_range_rates.csv');

% Get the number of epochs
% num_epochs = size(table_pseudo_ranges, 1) - 1; % Subtract 1 for the header
num_epochs = min(k, size(table_pseudo_ranges, 1) - 1);

% Initialize arrays to store results for all epochs
all_positions = zeros(num_epochs, 3);
all_clock_offsets = zeros(num_epochs, 1);
all_velocities = zeros(num_epochs, 3);
all_clock_drifts = zeros(num_epochs, 1);

% Initialize initial position at Earth center
r_ea_e = [0; 0; 0];
% Initialize receiver clock offset
receiver_clock_offset = 0;
% Initialize user velocity and clock drift
v_ea_e = [0; 0; 0];
receiver_clock_drift = 0;

% Initialize tolerance for convergence
tolerance = 0.1; % 10cm tolerance

for epoch = 1:num_epochs
    % 从CSV文件中读取卫星编号
    satellite_numbers = table_pseudo_ranges{1, 2:end};
    % Extract pseudo-ranges for the current epoch
    time = (epoch - 1) * 0.5; % Time starts from 0, 间隔为0.5
    pseudo_ranges = table_pseudo_ranges{1 + epoch, 2:end};
    pseudo_range_rates = table_pseudo_range_rates{1 + epoch, 2:end};
    
    iteration = 0;
    while true
        iteration = iteration + 1;
        % Predicted state vector
        x_prev = [r_ea_e; receiver_clock_offset];
        x_v_prev = [v_ea_e; receiver_clock_drift];

        % Calculate satellite positions and other parameters
        r_es_e = zeros(3, numel(satellite_numbers));
        v_es_e = zeros(3, numel(satellite_numbers));
        r_as = zeros(1, numel(satellite_numbers));
        v_as = zeros(1, numel(satellite_numbers));
        u_as_e = zeros(3, numel(satellite_numbers));
        delta_z = zeros(numel(satellite_numbers), 1);
        delta_z_v = zeros(numel(satellite_numbers), 1);
        H = zeros(numel(satellite_numbers), 4);

        for j = 1:numel(satellite_numbers)
            satellite_number = satellite_numbers(j);
            [r_ej_e, v_ej_e] = Satellite_position_and_velocity(time, satellite_number);
            r_es_e(:, j) = r_ej_e;
            v_es_e(:, j) = v_ej_e;
            
            C_I_e = eye(3);
            r_aj = norm(C_I_e * r_es_e(:, j) - r_ea_e);
            C_I_e = eye(3) - Omega_ie * r_aj / c;
            r_aj = norm(C_I_e * r_es_e(:, j) - r_ea_e);
            r_as(:, j) = r_aj;

            u_aj_e = (C_I_e * r_es_e(:, j) - r_ea_e) / r_aj;
            u_as_e(:, j) = u_aj_e;

            v_aj = u_aj_e' * (C_I_e * (v_ej_e' + Omega_ie * r_ej_e') - (v_ea_e + Omega_ie * r_ea_e));
            v_as(:, j) = v_aj;

            % Calculate measurements and measurement innovation
            delta_z(j) = pseudo_ranges(j) - (r_aj + receiver_clock_offset);
            delta_z_v(j) = pseudo_range_rates(j) - (v_aj + receiver_clock_drift);
            u_aj_e = u_as_e(:, j);
            H(j, :) = [-u_aj_e', 1];
        end

        % Calculate Moore-Penrose pseudo-inverse of H
        H_pseudo_inv = pinv(H);

        % Outlier Detection
        % 初始位置通常是一个估计值，会导致很大的残差。跳过第一个epoch以避免将初始估计的位置误认为是异常值。
        if epoch > 1
            % Calculate residuals
            residuals = (H * H_pseudo_inv - eye(numel(satellite_numbers))) * delta_z;
            % Calculate residuals covariance matrix
            Cv = (eye(numel(satellite_numbers)) - H * H_pseudo_inv) * (sigma_rho^2);
            % Check for outliers
            outliers = abs(residuals) > sqrt(diag(Cv)) * threshold;
            % Identify and remove outliers
            outlier_indices = find(outliers, 1);
            if ~isempty(outlier_indices)
                disp(['At time ' num2str(time) 's']);
                disp('Outlier Detected!');
                % Remove the measurement with the largest residual
                [~, max_residual_index] = max(abs(residuals));
                disp(['Remove outlier with the largest residual: Satellite ' num2str(satellite_numbers(max_residual_index))]);
                satellite_numbers(max_residual_index) = [];
                pseudo_ranges(max_residual_index) = [];
                continue;
            end
        end
        
        % Calculate the new state vector using least-squares
        x_new = x_prev + H_pseudo_inv * delta_z;
        % Extract position and receiver clock offset
        r_ea_e = x_new(1:3);
        receiver_clock_offset = x_new(4);

        x_v_new = x_v_prev + H_pseudo_inv * delta_z_v;
        v_ea_e = x_v_new(1:3);
        receiver_clock_drift = x_v_new(4);

        % Calculate the position difference between iterations
        position_difference = norm(r_ea_e - x_prev(1:3));
        % Check for convergence
        if position_difference < tolerance
            % disp(['经过', num2str(iteration), '次迭代']);
            break;
        end
        
    end

    % Store the results for this epoch
    all_positions(epoch, :) = r_ea_e';
    all_clock_offsets(epoch) = receiver_clock_offset;
    all_velocities(epoch, :) = v_ea_e';
    all_clock_drifts(epoch) = receiver_clock_drift;
end

% %%
% % Display the results
% disp('Time(s)   Latitude(deg)   Longitude(deg)   Height(m)   North(m/s)   East(m/s)   Down(m/s)');
% for epoch = 1:num_epochs
%     time = (epoch - 1) * 0.5;
%     position = all_positions(epoch, :);
%     velocity = all_velocities(epoch, :);
%     % Convert ECEF positions to latitude, longitude, and height
%     % Convert ECEF velocity to NED velocity
%     [L_b, lambda_b, h_b, v_eb_n] = pv_ECEF_to_NED(position', velocity');
%     % Convert latitude and longitude from radians to degrees
%     L_b_deg = rad_to_deg * L_b;
%     lambda_b_deg = rad_to_deg * lambda_b;
%     disp([num2str(time), '       ', num2str(L_b_deg), '        ', ...
%         num2str(lambda_b_deg), '       ', num2str(h_b), '       ', ...
%         num2str(v_eb_n(1)), '       ', num2str(v_eb_n(2)), '       ', num2str(v_eb_n(3))]);
% end
% 
% %%
% % Save results to a CSV file
% filename = 'results_GNSS_LSE.csv';
% % Open the CSV file for writing
% fileID = fopen(filename, 'w');
% % Check if the file was opened successfully
% if fileID == -1
%     error('Failed to open the file for writing');
% end
% 
% % Define the header as a cell array of strings
% header = {'Time(s)', 'Latitude(deg)', 'Longitude(deg)', 'Height(m)', 'North(m/s)', 'East(m/s)', 'Down(m/s)'};
% % Write the header to the CSV file
% fprintf(fileID, '%s,', header{1:end-1}); % Write all header elements except the last one
% fprintf(fileID, '%s\n', header{end});     % Write the last header element and a newline
% 
% % Write the data to the CSV file
% for epoch = 1:num_epochs
%     time = (epoch - 1) * 0.5;
%     position = all_positions(epoch, :);
%     velocity = all_velocities(epoch, :);
%     % Convert ECEF positions to latitude, longitude, and height
%     % Convert ECEF velocity to NED velocity
%     [L_b, lambda_b, h_b, v_eb_n] = pv_ECEF_to_NED(position', velocity');
%     % Convert latitude and longitude from radians to degrees
%     L_b_deg = rad_to_deg * L_b;
%     lambda_b_deg = rad_to_deg * lambda_b;
% 
%     data = [time, L_b_deg, lambda_b_deg, h_b, v_eb_n'];
% 
%     % Write the data to the CSV file
%     fprintf(fileID, '%f,', data(1:end-1)); % Write all data elements except the last one
%     fprintf(fileID, '%f\n', data(end));     % Write the last data element and a newline
% end
% 
% % Close the CSV file
% fclose(fileID);
% disp(['Results saved to "', filename, '"']);