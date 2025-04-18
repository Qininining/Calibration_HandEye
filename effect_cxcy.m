clc;
clear;
format longG;
%% 读取Excel文件中的数据
filename = 'Data.xlsx';    % Excel文件名
sheet = 1;                       % 工作表索引
try
    % 读取整列像素坐标数据
    points_pix = readmatrix(filename, 'Sheet', sheet, 'Range', 'A:B');
    % 读取整列世界坐标数据
    points_world = readmatrix(filename, 'Sheet', sheet, 'Range', 'E:G');
catch ME
    fprintf('读取Excel文件 "%s" 时出错: %s\n', filename, ME.message);
    fprintf('请确保文件存在，格式正确，并且包含所需的数据范围。\n');
    return; % 文件读取失败则退出
end
% 检查读取的数据是否为空
if isempty(points_pix) || isempty(points_world)
    fprintf('错误：从Excel读取的数据为空，请检查文件内容和范围设置。\n');
    return;
end
% 获取数据的行数
n = size(points_pix, 1);
if size(points_world, 1) ~= n
    fprintf('错误：像素坐标和世界坐标的点数不匹配 (%d vs %d)。\n', n, size(points_world, 1));
    return;
end


%% 验证dx dy的影响
dx_rand = randi([0, 1000]); % 随机生成dx，范围0-100,要求是integer
dy_rand = randi([0, 1000]); % 随机生成dy，范围0-100

fprintf('随机生成的dx: %d, dy: %d\n', dx_rand, dy_rand); % 修改格式化输出为整数

%% 转换像素坐标到相机坐标系下
% ================== 参数设置 ==================
% 标定参数（根据实际标定结果修改）
pixel_per_um_x = 820;      % X方向：1像素=0.82微米 (nm/px)
pixel_per_um_y = 820;      % Y方向：同上
dx = 800 + dx_rand;                   % 主点偏移x (px)
dy = 600 + dy_rand;                   % 主点偏移y (px)
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

%% 构建目标函数，求解最优参数 a, tx, ty, tz (使用 lsqnonlin)
% 将中间坐标转换为齐次坐标 (4xn)
points_cam_homogeneous = [points_cam, ones(n, 1)]'; % 齐次坐标

% 提取真实世界坐标 (nx3 -> 3xn)
points_world_real = points_world';

% 定义残差函数 (用于 lsqnonlin)
% params = [a, tx, ty, tz]
residual_func = @(params) calculate_residuals(params, points_cam_homogeneous, points_world_real);

% 设置初始猜测值
% 期望的旋转矩阵 cos(a) ≈ -0.9948, a ≈ +/- 3.0 rad
initial_guess = [-3.0, 0, 0, 0];  % 初始角度a (rad) 和 位移tx, ty, tz (nm)
fprintf('使用初始猜测值: a=%.4f rad, tx=%.1f, ty=%.1f, tz=%.1f\n', ...
        initial_guess(1), initial_guess(2), initial_guess(3), initial_guess(4));

% 设置优化选项
options = optimoptions('lsqnonlin', ...
    'Algorithm', 'levenberg-marquardt', ...
    'Display', 'iter-detailed', ... % 显示更详细的迭代信息
    'MaxFunctionEvaluations', 20000, ... % 增加评估次数上限
    'StepTolerance', 1e-9, ...      % 调整容差
    'FunctionTolerance', 1e-9);     % 调整容差
    % 'SpecifyObjectiveGradient',true, 'CheckGradients',false, ...
    % 'JacobianFcn', @(params) calculate_jacobian(params, points_cam_homogeneous) % (可选)

fprintf('开始使用 lsqnonlin 进行优化...\n');
% 求解优化问题
[params_optimized, resnorm, residual, exitflag, output] = lsqnonlin(residual_func, initial_guess, [], [], options);

% 检查优化是否成功
if exitflag <= 0
    fprintf('\n警告：优化可能未收敛或遇到问题。\n');
    disp(output); % 显示详细的输出信息
    % 即使优化失败，也尝试继续进行可视化以帮助诊断
else
    fprintf('\n优化成功完成。\n');
%     disp(output);
end

% 解析优化结果
a_optimized = params_optimized(1);
tx_optimized = params_optimized(2);
ty_optimized = params_optimized(3);
dz_optimized = params_optimized(4); % 对应参数 tz

% 显示优化结果
% fprintf('\nOptimized parameters:\n');
% fprintf('a  = %.8f (rad) / %.4f (deg)\n', a_optimized, rad2deg(a_optimized));
% fprintf('tx = %.6f (nm)\n', tx_optimized);
% fprintf('ty = %.6f (nm)\n', ty_optimized);
% fprintf('dz = %.6f (nm)\n', dz_optimized);
% fprintf('最终残差平方和 (resnorm): %.6e\n', resnorm);

% 使用优化后的参数构建最终的变换矩阵 T (从相机坐标系到世界坐标系)
cos_a = cos(a_optimized);
sin_a = sin(a_optimized);

% 注意 T 的结构，特别是 T(3,3) 和 T(3,4)
T_final = [cos_a,  -sin_a,   0,   tx_optimized;
           sin_a,   cos_a,   0,   ty_optimized;
           0,       0,      -1,   dz_optimized; % 这里的 -1 来自原始模型/期望结果
           0,       0,       0,   1];

disp('Final transformation matrix T (Camera Frame -> World Frame):');
disp(T_final);

%% 计算最终的估计世界坐标
points_world_estimated_final_homo = T_final * points_cam_homogeneous;
points_world_estimated_final = points_world_estimated_final_homo(1:3, :)'; % 转换回 nx3

% 计算最终的误差和均方根误差 (RMSE)
errors = points_world - points_world_estimated_final; % nx3 matrix, each row is [err_x, err_y, err_z]
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

T_M = T_final * M; % 计算最终的变换矩阵 T 和 M 的乘积
disp('Transformation matrix T * M:');
disp(T_M);

%% ================== 辅助函数 ==================
% 残差计算函数
function residuals = calculate_residuals(params, points_cam_homo, points_world_real)
    % params: [a, tx, ty, tz]
    a = params(1);
    tx = params(2);
    ty = params(3);
    tz = params(4); % 对应外部的 dz_optimized

    cos_a = cos(a);
    sin_a = sin(a);

    % 构建当前参数下的变换矩阵 T
    T = [cos_a,  -sin_a,   0,   tx;
         sin_a,   cos_a,   0,   ty;
         0,       0,      -1,   tz; % 这里的 -1 来自原始模型/期望结果
         0,       0,       0,   1];

    % 计算估计的世界坐标 (齐次，4xn)
    points_world_estimated_homo = T * points_cam_homo;

    % 提取估计的非齐次世界坐标 (3xn)
    points_world_estimated = points_world_estimated_homo(1:3, :);

    % 计算残差向量 (真实值 - 估计值)
    diff = points_world_real - points_world_estimated; % 3xn
    residuals = diff(:); % 返回列向量
end
