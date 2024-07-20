% Basic Kalman Filter Multiple Epochs

% load .csv data files
filename1 = 'Workshop2_GNSS_Pos_ECEF.csv';
pos_measurements = readmatrix(filename1);
% 总epoch数
num_epochs = size(pos_measurements, 1);

% Initialise the Kalman filter state vector estimate
% 初始估计位置 (x, y, z) 和速度 (vx, vy, vz) in ECEF system
x_ea_e = 2447019; 
y_ea_e = -5884199; 
z_ea_e = -284783; 
v_x_ea_e = 184;    
v_y_ea_e = 77;     
v_z_ea_e = 0;     
% 初始化状态向量 x
x_est = [x_ea_e; y_ea_e; z_ea_e; v_x_ea_e; v_y_ea_e; v_z_ea_e];
% Initialize arrays to store results for all epochs
x_ests = zeros(size(x_est, 1), num_epochs);

% Initialise the error covariance to give initial state uncertainties
% 初始不确定性的标准差
pos_uncertainty = 10;
vel_uncertainty = 5; 
% 初始化误差协方差矩阵 P
P_est = diag([pos_uncertainty^2, pos_uncertainty^2, pos_uncertainty^2, ...
                 vel_uncertainty^2, vel_uncertainty^2, vel_uncertainty^2]);

% (Kalman filter Step 1) Compute the transition matrix
t_s = 1;  % 给定传播间隔（秒）
I_3 = eye(3);
zero_3 = zeros(3);
% 使用分块矩阵创建状态转移矩阵 Phi
Phi = [I_3, t_s * I_3; zero_3, I_3];

% (Step 2) Compute the system noise covariance matrix
% 给定加速度功率谱密度（m^2/s^-3）
S_a_e = 5;% 计算系统噪声协方差矩阵 Q
Q = S_a_e * [1/3 * t_s^3 * I_3, 1/2 * t_s^2 * I_3; ...
             1/2 * t_s^2 * I_3, t_s * I_3];

for epoch = 1:num_epochs
    % (Step 3) Use the transition matrix to propagate the state estimates:
    x_predicted = Phi * x_est;
    
    % (Step 4) Then use transition matrix Phi & system noise covariance matrix Q
    % to propagate the error covariance matrix P:
    P_predicted = Phi * P_est * Phi' + Q;

    % (Step 5) Formulate the measurement matrix:
    H = [I_3, zero_3];
    
    % (Step 6) Compute the measurement noise covariance matrix 
    % assuming each component of the position measurement has an error standard deviation of 2.5 m
    sigma_x = 2.5;
    R = sigma_x^2 * I_3;

    % (Step 7) Compute the Kalman gain matrix
    % K = P_predicted * H' * inv(H * P_predicted * H' + R);
    K = P_predicted * H' / (H * P_predicted * H' + R);
    
    % (Step 8) Formulate the measurement innovation vector
    % Extract pos measurement 
    z = pos_measurements(epoch, 2:end)';
    z_innovation = z - x_predicted(1:3);
    
    % (Step 9) Update the state estimates
    x_est = x_predicted + K * z_innovation;
    
    % (Step 10) Update the error covariance matrix
    P_est = (eye(6) - K * H) * P_predicted;

    % Store the results for this epoch
    x_ests(:, epoch) = x_est;
    
end

% Display the results
disp('Time(s)   Latitude(deg)   Longitude(deg)   Height(m)   North(m/s)   East(m/s)   Down(m/s)');
for epoch = 1:num_epochs
    time = epoch - 1;
    x_est = x_ests(:, epoch);
    % Convert this Cartesian ECEF position solution to latitude, longitude and height
    [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(x_est(1:3), x_est(4:6));
    % 将纬度和经度从弧度转换为度
    L_b_deg = L_b * rad_to_deg;
    lambda_b_deg = lambda_b * rad_to_deg;
    disp([num2str(time), '       ', num2str(L_b_deg), '        ', ...
        num2str(lambda_b_deg), '       ', num2str(h_b), '       ', ...
        num2str(v_eb_n(1)), '       ', num2str(v_eb_n(2)), '       ', num2str(v_eb_n(3))]);
end


% Define the filename for the CSV file
filename = 'results_Task1b.csv';
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
