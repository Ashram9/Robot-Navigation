% First-epoch Positioning with Initialisation

Define_Constants

% load .csv data files
filename1 = 'Workshop1_Pseudo_ranges.csv';
table_pseudo_ranges = readtable(filename1);

% 定义移动电话的初始位置（纬度、经度和高度）
latitude_deg = -33.821075; % 纬度（度）
longitude_deg = 151.188496; % 经度（度）
height_m = 120; % 高度（米）
% 将纬度和经度转换为弧度
latitude_rad = latitude_deg * deg_to_rad;
longitude_rad = longitude_deg * deg_to_rad;
% 设置初始速度为零
initial_velocity = [0; 0; 0]; % 北、东、下方向的速度（米/秒）

% 调用pv_NED_to_ECEF函数进行坐标转换
[r_ea_e, ~] = pv_NED_to_ECEF(latitude_rad, longitude_rad, height_m, initial_velocity);
disp('手机初始ECEF坐标系下的位置：');
disp(r_ea_e);


% 从CSV文件中读取卫星编号
satellite_numbers = table_pseudo_ranges(1, 2:end);
% 初始化存储卫星位置的数组
r_es_e = zeros(3, numel(satellite_numbers));
% 设定时间为0s
time = 0;

% 计算time时刻下每颗卫星的ECEF坐标位置
for j = 1:numel(satellite_numbers)
    satellite_number = satellite_numbers.(j);
    [r_ej_e, ~] = Satellite_position_and_velocity(time, satellite_number);
    r_es_e(:, j) = r_ej_e;
end
% 输出计算得到的卫星位置
fprintf('Time = %ds, 卫星在ECEF坐标系下的位置：\n', time);
disp(r_es_e);


% 初始化存储估计距离的数组
r_as = zeros(1, numel(satellite_numbers));
% 初始化存储视线单位矢量的数组
u_as_e = zeros(3, numel(satellite_numbers));

for j = 1:numel(satellite_numbers)
    r_aj = norm(r_es_e(:, j) - r_ea_e);
    % 计算Sagnac效应补偿
    C_I_e = eye(3) - Omega_ie * r_aj / c;
    % 计算估计距离
    r_aj = norm(C_I_e * r_es_e(:, j) - r_ea_e);
    r_as(:, j) = r_aj;
    % 计算视线单位矢量
    u_aj_e = (C_I_e * r_es_e(:, j) - r_ea_e) / r_aj;
    u_as_e(:, j) = u_aj_e;
end
% 输出估计距离
disp('估计距离到各颗卫星的距离：');
disp(r_as);
% 输出视线单位矢量
disp('每颗卫星到用户位置的视线单位矢量：');
disp(u_as_e);


% predicted receiver clock offset
receiver_clock_offset = 0;
% predicted state vector
x = [r_ea_e; receiver_clock_offset];
% 初始化measurement innovation vector
delta_z = zeros(numel(satellite_numbers), 1);
% 初始化measurement matrix H
H = zeros(numel(satellite_numbers), 4);

% 从CSV文件中读取伪距
pseudo_ranges = table_pseudo_ranges(2+time, 2:end);
% 计算每颗卫星的测量量
for j = 1:numel(satellite_numbers)
    % 计算估计距离
    r_aj = r_as(j);
    % 计算测量创新
    delta_z(j) = pseudo_ranges.(j) - (r_aj + receiver_clock_offset);

    % 计算卫星到用户的视线单位矢量
    u_aj_e = u_as_e(:, j);
    % 构建测量矩阵 H
    H(j, :) = [-u_aj_e', 1];
end

disp('Predicted State Vector:');
disp(x);
disp('Measurement Innovation Vector:');
disp(delta_z);
disp('Measurement Matrix:');
disp(H);


% 计算H的伪逆矩阵（Moore-Penrose伪逆）
H_pseudo_inv = pinv(H);
% 计算最小二乘解 x
x_new = x + H_pseudo_inv * delta_z;

% 提取位置解和接收机钟偏差
receiver_position = x_new(1:3);
receiver_clock_offset = x_new(4);

% 输出位置解和接收机钟偏差
disp('Mobile Phone Position（ECEF坐标系下）：');
disp(receiver_position);
disp('Receiver Clock Offset（s）：');
disp(receiver_clock_offset);


% 使用pv_ECEF_to_NED函数将ECEF位置解转换为纬度、经度和高度
[L_b,lambda_b,h_b,~] = pv_ECEF_to_NED(receiver_position, [0; 0; 0]);

% 将纬度和经度从弧度转换为度
L_b_deg = L_b * rad_to_deg;
lambda_b_deg = lambda_b * rad_to_deg;

% 输出转换后的结果
disp('纬度（度）：');
disp(L_b_deg);
disp('经度（度）：');
disp(lambda_b_deg);
disp('高度（米）：');
disp(h_b);

