% First-epoch Positioning without Initialisation

Define_Constants

% load .csv data files
filename1 = 'Workshop1_Pseudo_ranges.csv';
table_pseudo_ranges = readtable(filename1);

% 定义移动电话的初始位置为地球中心(ECEF)
r_ea_e = [0; 0; 0];
% predicted receiver clock offset
receiver_clock_offset = 0;

% 设定时间为0min
time = 0;

% 从CSV文件中读取卫星编号
satellite_numbers = table_pseudo_ranges{1, 2:end};
% 从CSV文件中读取伪距
pseudo_ranges = table_pseudo_ranges{2+time, 2:end};


% 迭代参数
tolerance = 0.1; % 10cm的容忍度

iteration = 1;
while true
    % predicted state vector
    x_prev = [r_ea_e; receiver_clock_offset];

    % 初始化存储卫星位置的数组
    r_es_e = zeros(3, numel(satellite_numbers));
    % 初始化存储估计距离的数组
    r_as = zeros(1, numel(satellite_numbers));
    % 初始化存储视线单位矢量的数组
    u_as_e = zeros(3, numel(satellite_numbers));
    % 初始化measurement innovation vector
    delta_z = zeros(numel(satellite_numbers), 1);
    % 初始化measurement matrix H
    H = zeros(numel(satellite_numbers), 4);

    for j = 1:numel(satellite_numbers)
        satellite_number = satellite_numbers(j);
        % 计算time时刻下每颗卫星的ECEF坐标位置(记得将time从min转成s！！)
        [r_ej_e, ~] = Satellite_position_and_velocity(time*60, satellite_number);
        r_es_e(:, j) = r_ej_e;
    
        r_aj = norm(r_es_e(:, j) - r_ea_e);
        % 计算Sagnac效应补偿
        C_I_e = eye(3) - Omega_ie * r_aj / c;
        % 计算估计距离
        r_aj = norm(C_I_e * r_es_e(:, j) - r_ea_e);
        r_as(:, j) = r_aj;

        % 计算视线单位矢量
        u_aj_e = (C_I_e * r_es_e(:, j) - r_ea_e) / r_aj;
        u_as_e(:, j) = u_aj_e;

        % 计算测量创新
        delta_z(j) = pseudo_ranges(j) - (r_aj + receiver_clock_offset);
    
        % 计算卫星到用户的视线单位矢量
        u_aj_e = u_as_e(:, j);
        % 构建测量矩阵 H
        H(j, :) = [-u_aj_e', 1];
    end
    
    % 计算H的伪逆矩阵（Moore-Penrose伪逆）
    H_pseudo_inv = pinv(H);
    % 计算最小二乘解 x
    x_new = x_prev + H_pseudo_inv * delta_z;
    
    % 提取位置解和接收机钟偏差
    r_ea_e = x_new(1:3);
    receiver_clock_offset = x_new(4);
    
    % 计算位置解与上一次迭代的解之间的距离
    position_difference = norm(r_ea_e - x_prev(1:3));
    
    % 如果距离小于容忍度，停止迭代
    if position_difference < tolerance
        disp(['迭代收敛，迭代次数：' num2str(iteration)]);
        break;
    end
    
    iteration = iteration + 1;
end


% 输出位置解和接收机钟偏差
disp('Mobile Phone Position（ECEF坐标系下）：');
disp(r_ea_e);
disp('Receiver Clock Offset（s）：');
disp(receiver_clock_offset);


% 使用pv_ECEF_to_NED函数将ECEF位置解转换为纬度、经度和高度
[L_b,lambda_b,h_b,~] = pv_ECEF_to_NED(r_ea_e, [0; 0; 0]);

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
