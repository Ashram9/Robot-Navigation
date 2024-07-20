% Basic Kalman Filter State Uncertainties

% By taking the square roots of the 1st and 4th diagonal elements of the error covariance matrix,
% examine how the uncertainty of the x-direction position and velocity states vary over time.
% Note that they converge to constant values after a few seconds. 
% Repeat this exercise with the measurement update (steps 5 to 10) disabled and observe 
% how the uncertainties increase over time in the absence of new measurement information.

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

% Initialise the error covariance to give initial state uncertainties
% 初始不确定性的标准差
pos_uncertainty = 10;
vel_uncertainty = 5; 
% 初始化误差协方差矩阵 P
P_est = diag([pos_uncertainty^2, pos_uncertainty^2, pos_uncertainty^2, ...
                 vel_uncertainty^2, vel_uncertainty^2, vel_uncertainty^2]);
% 初始化用于记录不确定性随时间变化的数组
pos_uncertainties = zeros(1, num_epochs);
vel_uncertainties = zeros(1, num_epochs);

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
    % x_est = x_predicted;
    
    % (Step 10) Update the error covariance matrix
    P_est = (eye(6) - K * H) * P_predicted;
    % P_est = P_predicted;
    
    % log how the uncertainty of the x-direction position and velocity states vary over time
    pos_uncertainties(epoch) = sqrt(P_est(1, 1));
    vel_uncertainties(epoch) = sqrt(P_est(4, 4));
    
end

% 绘制x方向位置和速度状态的不确定性随时间的变化曲线
time_axis = 0:t_s:(num_epochs - t_s);
figure;
plot(time_axis, pos_uncertainties, 'b', 'LineWidth', 2);
hold on; % 在同一个窗口绘图
plot(time_axis, vel_uncertainties, 'r', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Uncertainty');
legend('Position Uncertainty', 'Velocity Uncertainty');
title('Uncertainty vs. Time');
grid on;  % 打开网格线

% The uncertianty of state (x-position & x-velocity) decrease and then converge to a constant over time.
% This is because the measurement information offset equals to the system noise.
% The uncertainty of position is smaller as the information from measurement is position.

% If we disable the measurement update, we can notice that the uncertainty increase dramatically.
% Besides, the rate of growth of position uncertianty is much  faster than velocity,
% as the position superimpose the uncertainty of itself as well as velocity.
