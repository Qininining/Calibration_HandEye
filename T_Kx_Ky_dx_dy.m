clc;
clear;
format longG;

%% ================== 符号表示 T*M ==================
% 定义符号变量
syms Rz Rx Ry tx ty tz sx sy dx dy Zf_sym real

% 构建符号旋转矩阵 R (ZYX 欧拉角)
cRz = cos(Rz); sRz = sin(Rz);
cRx = cos(Rx); sRx = sin(Rx);
cRy = cos(Ry); sRy = sin(Ry);
R_sym = [cRz*cRy, cRz*sRy*sRx - sRz*cRx, cRz*sRy*cRx + sRz*sRx;
         sRz*cRy, sRz*sRy*sRx + cRz*cRx, sRz*sRy*cRx - cRz*sRx;
         -sRy,   cRy*sRx,            cRy*cRx           ];

% 构建符号变换矩阵 T (相机 -> 世界)
T_sym = [ R_sym(1,1), R_sym(1,2), R_sym(1,3), tx;
          R_sym(2,1), R_sym(2,2), R_sym(2,3), ty;
          R_sym(3,1), R_sym(3,2), R_sym(3,3), tz;
          0,          0,          0,          1];

% 构建符号变换矩阵 M (像素 -> 相机)
M_sym = [sx,  0,  0, -dx*sx;
         0,  sy,  0, -dy*sy;
         0,  0,  1,  Zf_sym; % 使用符号 Zf
         0,  0,  0,  1];

% 计算符号组合变换矩阵 T_M
T_M_sym = T_sym * M_sym;

fprintf('符号表达式 T * M (像素 -> 世界):\n');
disp(T_M_sym);
fprintf('----------------------------------------\n\n');
% 清除符号变量，避免后续数值计算冲突
clear Rz Rx Ry tx ty tz sx sy dx dy Zf_sym cRz sRz cRx sRx cRy sRy R_sym T_sym M_sym T_M_sym

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

%% 定义相机内参的固定部分和初始猜测值
% ================== 固定参数 ==================
Zf = 10000000;      % 设定物距 (nm) - 注意：这里假设所有点都在同一深度, 保持固定

% ================== 初始猜测值 ==================
% T 矩阵参数的初始猜测 [Rz, Rx, Ry, tx, ty, tz]
initial_T_params = [3.1, 0.0, 0.0, 7.99e6, -1.25e7, 2.57e6];
% M 矩阵参数的初始猜测 [sx, sy, dx, dy]
initial_M_params = [820, 820, 400, 300]; % 使用之前的固定值作为初始猜测

% 组合所有优化参数的初始猜测
initial_guess = [initial_T_params, initial_M_params]; % 总共 10 个参数

fprintf('使用初始猜测值:\n');
fprintf('  T: Rz=%.4f, Rx=%.4f, Ry=%.4f rad, tx=%.1f, ty=%.1f, tz=%.1f\n', ...
        initial_guess(1), initial_guess(2), initial_guess(3), initial_guess(4), initial_guess(5), initial_guess(6));
fprintf('  M: sx=%.4f, sy=%.4f (nm/px), dx=%.1f, dy=%.1f (px)\n', ...
        initial_guess(7), initial_guess(8), initial_guess(9), initial_guess(10));


%% 构建目标函数，求解最优参数 (优化 T 和 M 的部分参数)
% 准备齐次像素坐标 (4xn)
points_pix_homo = [points_pix, zeros(n,1), ones(n,1)]'; % 输入 [u, v, 0, 1]

% 提取真实世界坐标 (3xn)
points_world_real = points_world';

% 定义残差函数 (用于 lsqnonlin)
% params = [Rz, Rx, Ry, tx, ty, tz, sx, sy, dx, dy] - 10 parameters
% 该函数计算 T*M*points_pix_homo 与 points_world_real 之间的差异
% 注意：Zf 作为固定参数传入
residual_func = @(params) calculate_residuals_TM_extended(params, points_pix_homo, points_world_real, Zf);

% 设置优化选项
options = optimoptions('lsqnonlin', ...
    'Algorithm', 'levenberg-marquardt', ...
    'Display', 'iter-detailed', ... % 显示更详细的迭代信息
    'MaxFunctionEvaluations', 30000, ... % 可能需要更多次迭代
    'StepTolerance', 1e-10, ...      % 可能需要更小的容差
    'FunctionTolerance', 1e-10);     % 可能需要更小的容差

fprintf('开始使用 lsqnonlin 进行优化 (优化 T 和 M 的参数)...\n');
% 求解优化问题
[params_optimized, resnorm, residual, exitflag, output] = lsqnonlin(residual_func, initial_guess, [], [], options);

% 检查优化是否成功
if exitflag <= 0
    fprintf('\n警告：优化可能未收敛或遇到问题。\n');
    disp(output); % 显示详细的输出信息
else
    fprintf('\n优化成功完成。\n');
end

% 解析优化结果
% T 的参数
Rz_optimized = params_optimized(1);
Rx_optimized = params_optimized(2);
Ry_optimized = params_optimized(3);
tx_optimized = params_optimized(4);
ty_optimized = params_optimized(5);
tz_optimized = params_optimized(6);
% M 的参数
sx_optimized = params_optimized(7);
sy_optimized = params_optimized(8);
dx_optimized = params_optimized(9);
dy_optimized = params_optimized(10);

% 显示优化结果
fprintf('\nOptimized parameters for T (Camera -> World):\n');
fprintf('Rz = %.8f (rad) / %.4f (deg)\n', Rz_optimized, rad2deg(Rz_optimized));
fprintf('Rx = %.8f (rad) / %.4f (deg)\n', Rx_optimized, rad2deg(Rx_optimized));
fprintf('Ry = %.8f (rad) / %.4f (deg)\n', Ry_optimized, rad2deg(Ry_optimized));
fprintf('tx = %.6f (nm)\n', tx_optimized);
fprintf('ty = %.6f (nm)\n', ty_optimized);
fprintf('tz = %.6f (nm)\n', tz_optimized);
fprintf('\nOptimized parameters for M (Pixel -> Camera):\n');
fprintf('sx = %.8f (nm/px)\n', sx_optimized);
fprintf('sy = %.8f (nm/px)\n', sy_optimized);
fprintf('dx = %.6f (px)\n', dx_optimized);
fprintf('dy = %.6f (px)\n', dy_optimized);
fprintf('Zf = %.1f (nm) (fixed)\n', Zf); % Zf 是固定的
fprintf('\n最终残差平方和 (resnorm): %.6e\n', resnorm);

% 使用优化后的参数构建最终的 T 矩阵
cRz = cos(Rz_optimized); sRz = sin(Rz_optimized);
cRx = cos(Rx_optimized); sRx = sin(Rx_optimized);
cRy = cos(Ry_optimized); sRy = sin(Ry_optimized);
R_final = [cRz*cRy, cRz*sRy*sRx - sRz*cRx, cRz*sRy*cRx + sRz*sRx;
           sRz*cRy, sRz*sRy*sRx + cRz*cRx, sRz*sRy*cRx - cRz*sRx;
           -sRy,   cRy*sRx,            cRy*cRx           ];
T_final = [ R_final(1,1), R_final(1,2), R_final(1,3), tx_optimized;
            R_final(2,1), R_final(2,2), R_final(2,3), ty_optimized;
            R_final(3,1), R_final(3,2), R_final(3,3), tz_optimized;
            0,            0,            0,            1];

% 使用优化后的参数构建最终的 M 矩阵
M_final = [sx_optimized,  0,  0, -dx_optimized*sx_optimized;
           0,  sy_optimized,  0, -dy_optimized*sy_optimized;
           0,  0,  1,  Zf;
           0,  0,  0,  1];

disp('Final transformation matrix T (Camera Frame -> World Frame):');
disp(T_final);
disp('Final transformation matrix M (Pixel Frame -> Camera Frame):');
disp(M_final);

% 计算最终的组合变换矩阵 T_M_final
T_M_final = T_final * M_final;
disp('Final combined transformation matrix T_M (Pixel Frame -> World Frame):');
disp(T_M_final);

%% 计算最终的估计世界坐标 (使用 T_M_final)
points_world_estimated_final_homo = T_M_final * points_pix_homo;
points_world_estimated_final = points_world_estimated_final_homo(1:3, :)'; % 转换回 nx3

% 计算最终的误差和均方根误差 (RMSE)
errors = points_world - points_world_estimated_final; % nx3 matrix, each row is [err_x, err_y, err_z]
rmse = sqrt(mean(sum(errors.^2, 2)));
fprintf('最终 RMSE: %.6f (nm)\n', rmse);

%% ================== 输出每个点的误差 ==================
fprintf('\n--- 每个点的误差详情 (单位: nm) ---\n');
fprintf('%-6s %-18s %-18s %-18s %-18s\n', '点号', '误差 X', '误差 Y', '误差 Z', '误差大小');
point_error_magnitudes = sqrt(sum(errors.^2, 2)); % Calculate magnitude for each point (nx1 vector)
for i = 1:n % 显示所有点的误差
    fprintf('%-6d %-18.6f %-18.6f %-18.6f %-18.6f\n', ...
            i, errors(i,1), errors(i,2), errors(i,3), point_error_magnitudes(i));
end
fprintf('----------------------------------------\n');

%% ================== 辅助函数 ==================
% 新的残差计算函数 (优化 T 和 M 的参数)
function residuals = calculate_residuals_TM_extended(params, points_pix_homo, points_world_real, Zf)
    % params: [Rz, Rx, Ry, tx, ty, tz, sx, sy, dx, dy] - 10 parameters
    % T 参数
    Rz = params(1); Rx = params(2); Ry = params(3);
    tx = params(4); ty = params(5); tz = params(6);
    % M 参数
    sx = params(7); sy = params(8); dx = params(9); dy = params(10);

    % 构建当前参数下的 T 矩阵
    cRz = cos(Rz); sRz = sin(Rz);
    cRx = cos(Rx); sRx = sin(Rx);
    cRy = cos(Ry); sRy = sin(Ry);
    R = [cRz*cRy, cRz*sRy*sRx - sRz*cRx, cRz*sRy*cRx + sRz*sRx;
         sRz*cRy, sRz*sRy*sRx + cRz*cRx, sRz*sRy*cRx - cRz*sRx;
         -sRy,   cRy*sRx,            cRy*cRx           ];
    T_current = [ R(1,1), R(1,2), R(1,3), tx;
                  R(2,1), R(2,2), R(2,3), ty;
                  R(3,1), R(3,2), R(3,3), tz;
                  0,      0,      0,      1];

    % 构建当前参数下的 M 矩阵 (使用固定的 Zf)
    M_current = [sx,  0,  0, -dx*sx;
                 0,  sy,  0, -dy*sy;
                 0,  0,  1,  Zf;
                 0,  0,  0,  1];

    % 计算当前的组合变换矩阵 T_M
    T_M_current = T_current * M_current;

    % 计算估计的世界坐标 (齐次，4xn)
    points_world_estimated_homo = T_M_current * points_pix_homo;

    % 提取估计的非齐次世界坐标 (3xn)
    points_world_estimated = points_world_estimated_homo(1:3, :);

    % 计算残差向量 (真实值 - 估计值)
    diff = points_world_real - points_world_estimated; % 3xn
    residuals = diff(:); % 返回列向量
end