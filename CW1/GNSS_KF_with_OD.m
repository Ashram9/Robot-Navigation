% GNSS Kalman Filter Multiple Epochs
addpath('Coursework 1 Data/');
addpath('Coursework 1 Software for Students to Use/');
%Define Constants
Define_Constants

% load .csv data files
pseudo_ranges = readmatrix('Pseudo_ranges.csv');
pseudo_range_rates = readmatrix('Pseudo_range_rates.csv');

% Get the number of epochs
num_epochs = size(pseudo_ranges, 1) - 1; % Subtract 1 for the header

% a)    
% Initialise state vector x using the solution obtained by LSE iteration method
disp('Initialising state by LSE...');
[r_ea_e, v_ea_e, receiver_clock_offset, receiver_clock_drift] = GNSS_LSE(1);
x_est = [r_ea_e; v_ea_e; receiver_clock_offset; receiver_clock_drift];
% Initialise error covariance matrix
% P_est = diag([10^2, 10^2, 10^2, 0.1^2, 0.1^2, 0.1^2, 100000^2, 200^2]);
P_est = diag([0, 0, 0, 0, 0, 0, 100000^2, 200^2]);

%%
% b) 
% (Kalman filter Step 1) Compute the transition matrix
t_s = 0.5;  % interval time (s)
I_3 = eye(3);
zero_3 = zeros(3);
zero_3_1 = zeros(3, 1);
zero_1_3 = zeros(1, 3);
% Construct transition matrix Phi
Phi = [I_3, t_s * I_3, zero_3_1, zero_3_1;
       zero_3, I_3, zero_3_1, zero_3_1;
       zero_1_3, zero_1_3, 1, t_s;
       zero_1_3, zero_1_3, 0, 1];

% c) 
% (Step 2) Compute the system noise covariance matrix
% acceleration PSD S_a_e = 5（m^2/s^-3)
% clock phase PSD S_cPhi_a = 0.01（m^2/s^-3)
% clock frequency PSD S_cf_a = 0.04（m^2/s^-1)
S_a_e = 5;
S_cPhi_a = 0.01;
S_cf_a = 0.04;
Q = [1/3 * S_a_e * t_s^3 * I_3, 1/2 * S_a_e * t_s^2 * I_3, zero_3_1, zero_3_1;
     1/2 * S_a_e * t_s^2 * I_3, S_a_e * t_s * I_3, zero_3_1, zero_3_1;
     zero_1_3, zero_1_3, S_cPhi_a * t_s + 1/3 * S_cf_a * t_s^3, 1/2 * S_cf_a * t_s^2;
     zero_1_3, zero_1_3, 1/2 * S_cf_a * t_s^2, S_cf_a * t_s];

% d) 
% (Step 3) Use the transition matrix to propagate the state estimates:
x_predicted = Phi * x_est;
r_ea_e = x_predicted(1:3);
v_ea_e = x_predicted(4:6);
receiver_clock_offset = x_predicted(7);
receiver_clock_drift = x_predicted(8);

% e) 
% (Step 4) Then use transition matrix Phi & system noise covariance matrix Q
% to propagate the error covariance matrix P:
P_predicted = Phi * P_est * Phi' + Q;

% Initialize arrays to store results for all epochs
x_ests = zeros(size(x_est, 1), num_epochs);

for epoch = 1: num_epochs
    time = (epoch - 1) * 0.5;
    % State Propogation
    x_predicted = Phi * x_est;
    r_ea_e = x_predicted(1:3);
    v_ea_e = x_predicted(4:6);
    receiver_clock_offset = x_predicted(7);
    receiver_clock_drift = x_predicted(8);
    P_predicted = Phi * P_est * Phi' + Q;

    % Measurement Update Preparation
    % Get indices of satelite
    satellite_numbers = pseudo_ranges(1, 2:end);
    pseudo_ranges_pt = pseudo_ranges(1+epoch, 2:end);
    pseudo_range_rates_pt = pseudo_range_rates(1+epoch, 2:end);

    while true
        H = zeros(2*numel(satellite_numbers), 8);
        z_innovation = zeros(2*numel(satellite_numbers), 1);

        for j = 1:numel(satellite_numbers)
            satellite_number = satellite_numbers(j);
            [r_ej_e, v_ej_e] = Satellite_position_and_velocity(time, satellite_number); 
            % f) 
            % Predict the ranges
            C_I_e = eye(3);
            r_aj = norm(C_I_e * r_ej_e' - r_ea_e);
            C_I_e = eye(3) - Omega_ie * r_aj / c;
            r_aj = norm(C_I_e * r_ej_e' - r_ea_e);
            % g)
            u_aj_e = (C_I_e * r_ej_e' - r_ea_e) / r_aj;
            % h) 
            % Predict the range rates
            v_aj = u_aj_e' * (C_I_e * (v_ej_e' + Omega_ie * r_ej_e') - (v_ea_e + Omega_ie * r_ea_e));
            % i) 
            % (Step 5) Compute the measurement matrix
            H(j, :) = [-u_aj_e', zero_1_3, 1, 0];
            H(j + numel(satellite_numbers), :) = [zero_1_3, -u_aj_e', 0, 1];
            % l) 
            % (Step 8) Formulate the measurement innovation vector
            z_innovation(j) = pseudo_ranges_pt(j) - (r_aj + receiver_clock_offset);
            z_innovation(j + numel(satellite_numbers)) = pseudo_range_rates_pt(j) - (v_aj + receiver_clock_drift);
        end
        
        % Outlier Detection
        sigma_rho = 5; % Measurement error standard deviation
        threshold = 6; % Outlier detection threshold
        % Calculate residuals
        H_check = H(1:numel(satellite_numbers), [1,2,3,7]);
        H_pseudo_inv = pinv(H_check);
        residuals = (H_check * H_pseudo_inv - eye(numel(satellite_numbers))) * z_innovation(1:numel(satellite_numbers));
        % Calculate residuals covariance matrix
        Cv = (eye(numel(satellite_numbers)) - H_check * H_pseudo_inv) * (sigma_rho^2);
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
            pseudo_ranges_pt(max_residual_index) = [];
            pseudo_range_rates_pt(max_residual_index) = [];
        else
            break;
        end

    end

    % j)
    % (Step 6) Compute the measurement noise covariance matrix 
    sigma_range = 10;
    sigma_range_rate = 0.05;
    R = [sigma_range^2 * eye(numel(satellite_numbers)), zeros(numel(satellite_numbers), numel(satellite_numbers));
         zeros(numel(satellite_numbers), numel(satellite_numbers)), sigma_range_rate^2 * eye(numel(satellite_numbers))];
    % h) 
    % (Step 7) Compute the Kalman gain matrix
    % K = P_predicted * H' * inv(H * P_predicted * H' + R);
    K = P_predicted * H' / (H * P_predicted * H' + R);

    % Measurement Update
    % m) 
    % (Step 9) Update the state estimates
    x_est = x_predicted + K * z_innovation;
    % k)
    % (Step 10) Update the error covariance matrix
    P_est = (eye(8) - K * H) * P_predicted;

    % Store the results for this epoch
    x_ests(:, epoch) = x_est;
    
end


% Save results to a CSV file
filename = 'results_GNSS_KF_with_OD.csv';
% Open the CSV file for writing
fileID = fopen(filename, 'w');
% Check if the file was opened successfully
if fileID == -1
    error('Failed to open the file for writing');
end

% Define the header as a cell array of strings
header = {'Time(s)', 'Latitude(deg)', 'Longitude(deg)', 'Height(m)', 'North(m/s)', 'East(m/s)', 'Down(m/s)'};
% Write the header to the CSV file
fprintf(fileID, '%s,', header{1:end-1}); % Write all header elements except the last one
fprintf(fileID, '%s\n', header{end});     % Write the last header element and a newline

% Write the data to the CSV file
for epoch = 1:num_epochs
    time = (epoch - 1) * 0.5;
    x_est = x_ests(:, epoch);
    % Convert this Cartesian ECEF position solution to latitude, longitude and height
    [L_b, lambda_b, h_b, v_eb_n] = pv_ECEF_to_NED(x_est(1:3), x_est(4:6));
    % Convert latitude and longitude from radians to degrees
    L_b_deg = L_b * rad_to_deg;
    lambda_b_deg = lambda_b * rad_to_deg;
    data = [time, L_b_deg, lambda_b_deg, h_b, v_eb_n'];
    
    % Write the data to the CSV file
    fprintf(fileID, '%f,', data(1:end-1)); % Write all data elements except the last one
    fprintf(fileID, '%f\n', data(end));     % Write the last data element and a newline
end

% Close the CSV file
fclose(fileID);
disp(['Results saved to "', filename, '"']);