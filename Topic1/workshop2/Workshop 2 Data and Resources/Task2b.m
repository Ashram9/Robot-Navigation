% GNSS Kalman Filter Multiple Epochs
%Define Constants
Define_Constants

% load .csv data files
filename1 = 'Workshop2_Pseudo_ranges.csv';
pseudo_ranges = readmatrix(filename1);
filename2 = 'Workshop2_Pseudo_range_rates.csv';
pseudo_range_rates = readmatrix(filename2);

% Get the number of epochs
num_epochs = size(pseudo_ranges, 1) - 1; % Subtract 1 for the header

% a)    
% 初始化状态向量 x & 误差协方差矩阵 P
[x_est, P_est] = Initialise_GNSS_KF;

% b) 
% (Kalman filter Step 1) Compute the transition matrix
t_s = 1;  % 给定传播间隔（秒）
I_3 = eye(3);
zero_3 = zeros(3);
zero_3_1 = zeros(3, 1);
zero_1_3 = zeros(1, 3);
% 使用分块矩阵创建状态转移矩阵 Phi
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
    time = epoch - 1;
    % State Propogation
    x_predicted = Phi * x_est;
    r_ea_e = x_predicted(1:3);
    v_ea_e = x_predicted(4:6);
    receiver_clock_offset = x_predicted(7);
    receiver_clock_drift = x_predicted(8);
    P_predicted = Phi * P_est * Phi' + Q;

    % Measurement Update Preparation
    % 从CSV文件中读取卫星编号
    satellite_numbers = pseudo_ranges(1, 2:end);
    pseudo_ranges_pt = pseudo_ranges(1+epoch, 2:end);
    pseudo_range_rates_pt = pseudo_range_rates(1+epoch, 2:end);
    
    H = zeros(numel(satellite_numbers), 8);
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
    
    % j)
    % (Step 6) Compute the measurement noise covariance matrix 
    sigma_range = 10;
    sigma_range_rate = 0.05;
    I_10 = eye(10);
    zero_10 = zeros(10);
    R = [sigma_range^2 * I_10, zero_10;
         zero_10, sigma_range_rate^2 * I_10];
    
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

% o)
% Display the results
disp('Time(s)   Latitude(deg)   Longitude(deg)   Height(m)   North(m/s)   East(m/s)   Down(m/s)');
for epoch = 1:num_epochs
    time = epoch - 1;
    x_est = x_ests(:, epoch);
    position = x_est(1:3);
    velocity = x_est(4:6);
    % Convert ECEF positions to latitude, longitude, and height
    % Convert ECEF velocity to NED velocity
    [L_b, lambda_b, h_b, v_eb_n] = pv_ECEF_to_NED(position, velocity);
    % Convert latitude and longitude from radians to degrees
    L_b_deg = rad_to_deg * L_b;
    lambda_b_deg = rad_to_deg * lambda_b;
    disp([num2str(time), '       ', num2str(L_b_deg), '        ', ...
        num2str(lambda_b_deg), '       ', num2str(h_b), '       ', ...
        num2str(v_eb_n(1)), '       ', num2str(v_eb_n(2)), '       ', num2str(v_eb_n(3))]);
end

% Save results to a CSV file
filename = 'results_Task2b.csv';
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
    time = epoch - 1;
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