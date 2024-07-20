% Gyroscopic & Magnetic Integration

addpath('Coursework 1 Data/');
addpath('Coursework 1 Software for Students to Use/');
% Define_Constants
Define_Constants

% load .csv data file
data_DR = readmatrix('Dead_reckoning.csv');
times = data_DR(:, 1); % Time in seconds
angular_rate = data_DR(:, 6); % Gyroscopic angular rate
heading_M = data_DR(:, 7) * deg_to_rad; % Magnetic Heading

% store Gyroscopic Heading
heading_G = zeros(size(heading_M));
%% 
% Step1: calculate Gyroscopic Heading
% initialise
heading_G(1) = heading_M(1); 

for k = 2:length(times)
    heading_G(k) = heading_G(k-1) + angular_rate(k) * (times(k) - times(k-1));
end

%%
% Step2: Gyro-Magnetometer Integration
% store G&M integrated Heading
heading_G_M = zeros(size(heading_M));
heading_G_M(1) = heading_M(1);
% initialise
x_est = zeros(2, 1);
x_est(1) = heading_G(1) - heading_M(1);
% gyro bias is 1 degree per second
x_est(2) = 1 * deg_to_rad;
P_est = zeros(2);

t_s = 0.5;

for k = 2:length(times)
    
    % 预测步骤...
    % (Kalman filter Step 1) Compute the transition matrix
    Phi = [1, t_s;
           0, 1];

    % (Step 2) Compute the system noise covariance matrix
    % Q = zeros(2);
    S_rg = 0;
    S_bgd = 0;
    Q = [S_rg * t_s + 1/3 * S_bgd * t_s^3, 1/2 * S_bgd * t_s^2;
         1/2 * S_bgd * t_s^2, S_bgd * t_s];

    % (Step 3) Propagate the state estimates:
    x_predicted = Phi * x_est;
    % (Step 4) Propagate the error covariance matrix:
    P_predicted = Phi * P_est * Phi' + Q;
    
    % 更新步骤...
    % (Step 5) Compute the measurement matrix
    H = [-1, 0];
    % (Step 6) Compute the measurement noise covariance matrix
    % magnetix heading noise variance
    sigma_M = 4 * deg_to_rad;
    R = sigma_M^2;

    % (Step 7) Compute the Kalman gain matrix
    K = P_predicted * H' / (H * P_predicted * H' + R);

    % (Step 8) Formulate the measurement innovation vector
    % 测量量z = GNSS数据 - DR数据
    z = heading_M(k) - heading_G(k);
    z_innovation = z - H * x_predicted;
    
    % (Step 9) Update the state estimates
    x_est = x_predicted + K * z_innovation;
    % (Step 10) Update the error covariance matrix
    P_est = (eye(2) - K * H) * P_predicted;
    
    % Correct the DR solution
    heading_G_M(k) =  heading_G(k) - x_est(1);
end

%%
% Define the filename for the CSV file
filename = 'headings.csv';
% Open the CSV file for writing
fileID = fopen(filename, 'w');
% Check if the file was opened successfully
if fileID == -1
    error('Failed to open the file for writing');
end

% Define the header as a cell array of strings
header = {'Time(s)', 'Gyro Heading(deg)', 'Mag Heading(deg)', 'Integrated Heading(deg)'};
% Write the header to the CSV file
fprintf(fileID, '%s,', header{1:end-1}); % Write all header elements except the last one
fprintf(fileID, '%s\n', header{end});     % Write the last header element and a newline

% Write the data to the CSV file
for k = 1:length(times)
    time = times(k);
    data = [time, heading_G(k) * rad_to_deg, heading_M(k) * rad_to_deg, heading_G_M(k) * rad_to_deg];
    
    % Write the data to the CSV file
    fprintf(fileID, '%f,', data(1:end-1)); % Write all data elements except the last one
    fprintf(fileID, '%f\n', data(end));     % Write the last data element and a newline
end

% Close the CSV file
fclose(fileID);
disp(['Results saved to "', filename, '"']);