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
% dx_rand = randi([0, 1000]); % 随机生成dx，范围0-100,要求是integer
% dy_rand = randi([0, 1000]); % 随机生成dy，范围0-100
% 
% fprintf('随机生成的dx: %d, dy: %d\n', dx_rand, dy_rand); % 修改格式化输出为整数

%% 转换像素坐标到相机坐标系下
% ================== 参数设置 ==================
% 标定参数（根据实际标定结果修改）
pixel_per_um_x = 820;      % X方向：1像素=0.82微米 (nm/px)
pixel_per_um_y = 820;      % Y方向：同上
% dx = 800 + dx_rand;                   % 主点偏移x (px)
% dy = 600 + dy_rand;                   % 主点偏移y (px)
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

%% 构建目标函数，求解最优参数 a, tx, ty, tz (使用 lsqnonlin)
% 将中间坐标转换为齐次坐标 (4xn)
points_cam_homogeneous = [points_cam, ones(n, 1)]'; % 齐次坐标

% 提取真实世界坐标 (nx3 -> 3xn)
points_world_real = points_world';

% 定义残差函数 (用于 lsqnonlin)
% params = [Rz, Rx, Ry, tx, ty, tz] - ZYX Euler angles and translation
residual_func = @(params) calculate_residuals(params, points_cam_homogeneous, points_world_real);

% 设置初始猜测值
% 期望的旋转矩阵 cos(Rz) ≈ -0.9989, Rz ≈ +/- 3.1 rad. R(3,3) ≈ -1 suggests potential issues or different conventions.
% Let's use the translation from the expected result and angles close to 180 deg rotation around Z.
% Original guess: initial_guess = [-3.0, -3, 0, 0, 0, 7566225]; % Corresponded to [a, b, r, tx, ty, tz]
% Revised guess based on expected result and ZYX order [Rz, Rx, Ry, tx, ty, tz]:
initial_guess = [-3.1, 0.0, 0.0, 7.99e6, -1.25e7, 7.57e6]; % 使用更接近期望值的初始猜测 (Rz≈-177deg, Rx≈0, Ry≈0, tx, ty, tz from expected)
fprintf('使用修正后的初始猜测值: Rz=%.4f, Rx=%.4f, Ry=%.4f rad, tx=%.1f, ty=%.1f, tz=%.1f\n', ...
        initial_guess(1), initial_guess(2), initial_guess(3), initial_guess(4), initial_guess(5), initial_guess(6));

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
Rz_optimized = params_optimized(1); % Rotation around Z
Rx_optimized = params_optimized(2); % Rotation around X
Ry_optimized = params_optimized(3); % Rotation around Y
tx_optimized = params_optimized(4);
ty_optimized = params_optimized(5);
dz_optimized = params_optimized(6); % 对应参数 tz

% 显示优化结果
fprintf('\nOptimized parameters:\n');
fprintf('Rz = %.8f (rad) / %.4f (deg)\n', Rz_optimized, rad2deg(Rz_optimized));
fprintf('Rx = %.8f (rad) / %.4f (deg)\n', Rx_optimized, rad2deg(Rx_optimized));
fprintf('Ry = %.8f (rad) / %.4f (deg)\n', Ry_optimized, rad2deg(Ry_optimized));
fprintf('tx = %.6f (nm)\n', tx_optimized);
fprintf('ty = %.6f (nm)\n', ty_optimized);
fprintf('dz = %.6f (nm)\n', dz_optimized);
fprintf('最终残差平方和 (resnorm): %.6e\n', resnorm);

% 使用优化后的参数构建最终的变换矩阵 T (从相机坐标系到世界坐标系)
% Calculate sines and cosines for optimized angles
cRz = cos(Rz_optimized); sRz = sin(Rz_optimized);
cRx = cos(Rx_optimized); sRx = sin(Rx_optimized);
cRy = cos(Ry_optimized); sRy = sin(Ry_optimized);

% Rotation matrix R from ZYX Euler angles: R = Rz(Rz) * Ry(Ry) * Rx(Rx)
% Note: The variable names 'a', 'b', 'r' in the original formula correspond to Rz, Rx, Ry here.
% R = [ca*cr, ca*sr*sb - sa*cb, ca*sr*cb + sa*sb;
%      sa*cr, sa*sr*sb + ca*cb, sa*sr*cb - ca*sb;
%      -sr,   cr*sb,            cr*cb           ];
% Translating to new variables (a=Rz, b=Rx, r=Ry):
R_final = [cRz*cRy, cRz*sRy*sRx - sRz*cRx, cRz*sRy*cRx + sRz*sRx;
           sRz*cRy, sRz*sRy*sRx + cRz*cRx, sRz*sRy*cRx - cRz*sRx;
           -sRy,   cRy*sRx,            cRy*cRx           ];


% 构建最终的变换矩阵 T_final
T_final = [ R_final(1,1), R_final(1,2), R_final(1,3), tx_optimized;
            R_final(2,1), R_final(2,2), R_final(2,3), ty_optimized;
            R_final(3,1), R_final(3,2), R_final(3,3), dz_optimized;
            0,            0,            0,            1];

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
for i = 1%n
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
    % params: [Rz, Rx, Ry, tx, ty, tz] - ZYX Euler angles and translation
    Rz = params(1); % Rotation around Z
    Rx = params(2); % Rotation around X
    Ry = params(3); % Rotation around Y
    tx = params(4);
    ty = params(5);
    tz = params(6); % 对应外部的 dz_optimized

    % Calculate sines and cosines
    cRz = cos(Rz); sRz = sin(Rz);
    cRx = cos(Rx); sRx = sin(Rx);
    cRy = cos(Ry); sRy = sin(Ry);

    % Rotation matrix R from ZYX Euler angles: R = Rz(Rz) * Ry(Ry) * Rx(Rx)
    % Note: The variable names 'a', 'b', 'r' in the original formula correspond to Rz, Rx, Ry here.
    % R = [ca*cr, ca*sr*sb - sa*cb, ca*sr*cb + sa*sb;
    %      sa*cr, sa*sr*sb + ca*cb, sa*sr*cb - ca*sb;
    %      -sr,   cr*sb,            cr*cb           ];
    % Translating to new variables (a=Rz, b=Rx, r=Ry):
    R = [cRz*cRy, cRz*sRy*sRx - sRz*cRx, cRz*sRy*cRx + sRz*sRx;
         sRz*cRy, sRz*sRy*sRx + cRz*cRx, sRz*sRy*cRx - cRz*sRx;
         -sRy,   cRy*sRx,            cRy*cRx           ];


    % 构建当前参数下的变换矩阵 T
    T = [ R(1,1), R(1,2), R(1,3), tx;
          R(2,1), R(2,2), R(2,3), ty;
          R(3,1), R(3,2), R(3,3), tz;
          0,      0,      0,      1];

    % 计算估计的世界坐标 (齐次，4xn)
    points_world_estimated_homo = T * points_cam_homo;

    % 提取估计的非齐次世界坐标 (3xn)
    points_world_estimated = points_world_estimated_homo(1:3, :);

    % 计算残差向量 (真实值 - 估计值)
    diff = points_world_real - points_world_estimated; % 3xn
    residuals = diff(:); % 返回列向量
end
