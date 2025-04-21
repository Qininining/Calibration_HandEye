clc
clear
format longG;

% 读取Excel文件中的数据
filename = 'Data.xlsx';    % Excel文件名
sheet = 1;                       % 工作表索引

% 读取整列像素坐标数据
points_pix = readmatrix(filename, 'Sheet', sheet, 'Range', 'A:B');

% 获取数据的行数
n = height(points_pix);

% ================== 参数设置 ==================
% 标定参数（根据实际标定结果修改）
pixel_per_um_x = 820;      % X方向：1像素=0.82微米 (nm/px)
pixel_per_um_y = 820;      % Y方向：同上
dx = 800;                   % 主点偏移x (px)
dy = 600;                   % 主点偏移y (px)
Zf = 10000000;                 % 设定物距 (nm)

% ================== 构建转换矩阵M ==================
M = [pixel_per_um_x, 0,          -dx*pixel_per_um_x;
     0,               pixel_per_um_y, -dy*pixel_per_um_y;
     0,               0,           Zf];

% ================== 坐标转换 ==================
% 添加齐次坐标分量1
points_homo = [points_pix, ones(n, 1)];

% 矩阵乘法转换
points_cam = points_homo * M'; 

% % ================== 结果显示 ==================
% disp('原始像素坐标：');
% disp(points_pix);
% disp('转换后的相机坐标（单位：μm）：');
% disp(points_cam(:,1:3)); % 显示XYZ三维坐标

% 读取整列世界坐标数据
points_world = readmatrix(filename, 'Sheet', sheet, 'Range', 'E:G');



%% 
A = [points_cam(:,1:2), ones(size(points_cam(:,1:2),1),1)];

%% 3. 分离目标向量分量
b_x = points_world(:,1);
b_y = points_world(:,2);
b_z = points_world(:,3);


%% 5. 求解最小二乘解
m_x = A \ b_x;
m_y = A \ b_y;
m_z = A \ b_z;

%% 6. 构建变换矩阵
T = [
    m_x(1), m_x(2), 0,  m_x(3);
    m_y(1), m_y(2), 0,  m_y(3);
    0,      0,      -1, m_z(3) + Zf;
    0,      0,      0,  1;
];

% 
T_rounded = round(T, 6);
disp('计算得到的变换矩阵:');
disp(T_rounded);

% %% 7. 验证结果
% transformed_points = (T * [points_cam, ones(n,1)]')';
% transformed_points = transformed_points(:,1:3);
% error = norm(transformed_points - points_world, 'fro');
% fprintf('变换误差(Frobenius范数): %.6f nm\n', error);
% 
% %% 计算是否正交
% seita1 = asin(abs(m_x(2)/m_x(1))) / pi *180;
% seita2 = asin(abs(m_y(1)/m_y(2))) / pi *180;
% err = abs(seita1-seita2);
% 
% %% 正交 
% T = [
%     (m_x(1) + m_y(2))/2,         -(abs(m_x(2)) + m_y(1))/2,      0,     m_x(3);
%     (abs(m_x(2)) + m_y(1))/2,    (m_x(1) + m_y(2))/2,            0,     m_y(3);
%     0,                           0,                             -1,     m_z(3) + Zf;
%     0,                           0,                             0,      1;
% ];
% T_rounded = round(T, 6);
% disp('正交化后的变换矩阵:');
% disp(T_rounded);
% 
% %% 用角度表示
% a = asin(( (abs(m_x(2)) + m_y(1))/2)/((m_x(1) + m_y(2))/2)) + pi;
% aa = a/pi * 180 - 180;
% disp('角度偏差:(°)');
% disp(aa);
% % a = a
% % a = seita1 / 180 * pi + pi
% % a = seita2 / 180 * pi + pi
% 
% T = [
%     cos(a),         -sin(a),      0,     m_x(3);
%     sin(a),    cos(a),            0,     m_y(3);
%     0,                           0,                             -1,     m_z(3) + Zf;
%     0,                           0,                             0,      1;
% ];
% 
% T_rounded = round(T, 6);
% disp('用角度表示后的变换矩阵:');
% disp(T_rounded);

%% 单位化

%% 7. 验证结果
transformed_points = (T * [points_cam, ones(n,1)]')';
transformed_points = transformed_points(:,1:3);
error = norm(transformed_points - points_world, 'fro');
fprintf('变换误差(Frobenius范数): %.6f nm\n', error);



%% (可选) 计算最终的估计世界坐标并比较
% points_world_estimated_final_homo = T * points_cam_homogeneous;
points_world_estimated_final = transformed_points; % 转换回 nx3

% 计算最终的均方根误差 (RMSE)
errors = points_world - points_world_estimated_final;
rmse = sqrt(mean(sum(errors.^2, 2)));
fprintf('最终 RMSE: %.6f (nm)\n', rmse);
%% ================== 输出每个点的误差 ==================
fprintf('\n--- 每个点的误差详情 (单位: nm) ---\n');
fprintf('%-6s %-18s %-18s %-18s %-18s\n', '点号', '误差 X', '误差 Y', '误差 Z', '误差大小');
point_error_magnitudes = sqrt(sum(errors.^2, 2)); % Calculate magnitude for each point (nx1 vector)
for i = 1:n
    fprintf('%-6d %-18.6f %-18.6f %-18.6f %-18.6f\n', ...
            i, errors(i,1), errors(i,2), errors(i,3), point_error_magnitudes(i));
end
fprintf('----------------------------------------\n');


%% 计算T*M
% 将M转换为4*4
M = [pixel_per_um_x,    0,              0,  -dx*pixel_per_um_x;
     0,                 pixel_per_um_y, 0,  -dy*pixel_per_um_y;
     0,                 0,              1,  Zf;
     0,                 0,              0,  1]; % 添加齐次坐标分量

T_M = T * M; % 计算最终的变换矩阵 T 和 M 的乘积
disp('Transformation matrix T * M:');
disp(T_M);
