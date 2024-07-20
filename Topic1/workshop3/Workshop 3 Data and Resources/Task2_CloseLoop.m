% Car DR/GNSS Integration (Close-Loop)
% Define_Constants
Define_Constants

% 加载GNSS数据
data_GNSS = readmatrix('Workshop3_GNSS_Pos_Vel_NED.csv');
times = data_GNSS(:, 1);
latitudes_rad_GNSS = data_GNSS(:, 2) * deg_to_rad;
longitudes_rad_GNSS = data_GNSS(:, 3) * deg_to_rad;
height_GNSS = data_GNSS(:, 4);
v_N_GNSS = data_GNSS(:, 5);
v_E_GNSS = data_GNSS(:, 6);
% v_D_GNSS = data(:, 7);

% 加载DR数据（从Task1中计算得到）
data_DR = readmatrix('results_Task1.csv');
latitudes_rad_DR = data_DR(:, 2) * deg_to_rad;
longitudes_rad_DR = data_DR(:, 3) * deg_to_rad;
v_N_DR = data_DR(:, 4);
v_E_DR = data_DR(:, 5);

% load .csv data file
data = readmatrix('Workshop3_Speed_Heading.csv');
speed = data(:, 2); % Forward speed in m/s
heading_deg = data(:, 3); % Heading in degrees
heading_rad = heading_deg * deg_to_rad; % Heading in radians
% store average velocities
avg_velocities = zeros(length(times), 2);

% 初始状态向量和误差协方差矩阵
x_est = zeros(4, 1);
[RN, RE] = Radii_of_curvature(latitudes_rad_GNSS(1)); % Compute meridian and transverse radii of curvature
P_est = diag([0.1^2, 0.1^2, 10^2/(RN + height_GNSS(1))^2, 10^2/((RE + height_GNSS(1))^2 * cos(latitudes_rad_GNSS(1))^2)]);

t_s = 0.5;

% 初始化0时刻corrected DR solution 为 GNSS solution
latitudes_rad_DR(1) = latitudes_rad_GNSS(1);
longitudes_rad_DR(1) = longitudes_rad_GNSS(1);
v_N_DR(1) = v_N_GNSS(1);
v_E_DR(1) = v_E_GNSS(1);

avg_velocities(1, :) = [v_N_DR(1), v_E_DR(1)];

% 卡尔曼滤波器的预测和更新步骤
for k = 2:length(times)
    %% 
    % Close-Loop (recalculate DR)
    heading_cur = heading_rad(k);
    heading_prev = heading_rad(k-1);

    % Resolving Vector
    M_NE = 1/2 * [cos(heading_cur) + cos(heading_prev), sin(heading_cur) + sin(heading_prev)];
    avg_velocities(k, :) = M_NE * speed(k);
    v_N_DR(k) = 1.7 * avg_velocities(k, 1) - 0.7 * v_N_DR(k-1);
    v_E_DR(k) = 1.7 * avg_velocities(k, 2) - 0.7 * v_E_DR(k-1);

    [RN, RE] = Radii_of_curvature(latitudes_rad_DR(k-1)); % Compute meridian and transverse radii of curvature

    velocity_N = avg_velocities(k, 1);
    velocity_E = avg_velocities(k, 2);

    % Calculate latitude change using the formula
    delta_latitude = velocity_N * t_s / (RN + height_GNSS(k));
    % Update latitude
    latitudes_rad_DR(k) = latitudes_rad_DR(k-1) + delta_latitude;
    % Calculate longitude change using the formula
    delta_longitude = velocity_E * t_s / ((RE + height_GNSS(k)) * cos(latitudes_rad_DR(k)));
    % Update longitude
    longitudes_rad_DR(k) = longitudes_rad_DR(k-1) + delta_longitude;
    % 注意这里要清零！！因为已经校正了
    x_est = zeros(4, 1);
    %%
    %计算地球相关半径(根据前一时刻的latitudes_rad_GNSS)
    [RN, RE] = Radii_of_curvature(latitudes_rad_GNSS(k-1));
    
    % 预测步骤...
    % (Kalman filter Step 1) Compute the transition matrix
    Phi = [1, 0, 0, 0;
           0, 1, 0, 0;
           t_s / (RN + height_GNSS(k-1)), 0, 1, 0;
           0, t_s / ((RE + height_GNSS(k-1)) * cos(latitudes_rad_GNSS(k-1))), 0, 1];

    % (Step 2) Compute the system noise covariance matrix
    % DR velocity error PSD S_DR = 0.2（m^2/s^-3)
    S_DR = 0.2;
    Q = [S_DR * t_s, 0,  1/2 * S_DR * t_s^2 / (RN + height_GNSS(k-1)), 0;
         0, S_DR * t_s, 0, 1/2 * S_DR * t_s^2 / ((RE + height_GNSS(k-1)) * cos(latitudes_rad_GNSS(k-1)));
         1/2 * S_DR * t_s^2 / (RN + height_GNSS(k-1)), 0, 1/3 * S_DR * t_s^3 / (RN + height_GNSS(k-1))^2, 0;
         0, 1/2 * S_DR * t_s^2 / ((RE + height_GNSS(k-1)) * cos(latitudes_rad_GNSS(k-1))), 0, 1/3 * S_DR * t_s^3 / ((RE + height_GNSS(k-1))^2 * cos(latitudes_rad_GNSS(k-1))^2)];

    % (Step 3) Propagate the state estimates:
    x_predicted = Phi * x_est;
    % (Step 4) Propagate the error covariance matrix:
    P_predicted = Phi * P_est * Phi' + Q;
    
    % 更新步骤...
    % (Step 5) Compute the measurement matrix
    H = [0, 0, -1, 0;
         0, 0, 0, -1;
         -1, 0, 0, 0;
         0, -1, 0, 0];
    % (Step 6) Compute the measurement noise covariance matrix 
    sigma_Gr = 5;
    sigma_Gv = 0.02;
    R = diag([sigma_Gr^2 / (RN + height_GNSS(k))^2, ...
              sigma_Gr^2 / ((RE + height_GNSS(k))^2 * cos(latitudes_rad_GNSS(k))^2), ...
              sigma_Gv^2, sigma_Gv^2]);

    % (Step 7) Compute the Kalman gain matrix
    K = P_predicted * H' / (H * P_predicted * H' + R);

    % (Step 8) Formulate the measurement innovation vector
    % 测量量z = GNSS数据 - DR数据
    z = [latitudes_rad_GNSS(k) - latitudes_rad_DR(k);
         longitudes_rad_GNSS(k) - longitudes_rad_DR(k);
         v_N_GNSS(k) - v_N_DR(k);
         v_E_GNSS(k) - v_E_DR(k)];
    z_innovation = z - H * x_predicted;
    
    % (Step 9) Update the state estimates
    x_est = x_predicted + K * z_innovation;
    % (Step 10) Update the error covariance matrix
    P_est = (eye(4) - K * H) * P_predicted;
    
    % Correct the DR solution
    latitudes_rad_DR(k) = latitudes_rad_DR(k) - x_est(3);
    longitudes_rad_DR(k) = longitudes_rad_DR(k) - x_est(4);
    v_N_DR(k) = v_N_DR(k) - x_est(1);
    v_E_DR(k) = v_E_DR(k) - x_est(2);
    

end

% Use the Kalman filter estimates to correct the DR solution at each epoch, k

%%
% Define the filename for the CSV file
filename = 'results_Task2_CloseLoop.csv';
% Open the CSV file for writing
fileID = fopen(filename, 'w');
% Check if the file was opened successfully
if fileID == -1
    error('Failed to open the file for writing');
end

% Define the header as a cell array of strings
header = {'Time(s)', 'Latitude(deg)', 'Longitude(deg)', 'North(m/s)', 'East(m/s)'};
% Write the header to the CSV file
fprintf(fileID, '%s,', header{1:end-1}); % Write all header elements except the last one
fprintf(fileID, '%s\n', header{end});     % Write the last header element and a newline

% Write the data to the CSV file
for k = 1:length(times)
    time = times(k);
    data = [time, latitudes_rad_DR(k) * rad_to_deg, longitudes_rad_DR(k) * rad_to_deg, v_N_DR(k), v_E_DR(k)];
    
    % Write the data to the CSV file
    fprintf(fileID, '%f,', data(1:end-1)); % Write all data elements except the last one
    fprintf(fileID, '%f\n', data(end));     % Write the last data element and a newline
end

% Close the CSV file
fclose(fileID);
disp(['Results saved to "', filename, '"']);


