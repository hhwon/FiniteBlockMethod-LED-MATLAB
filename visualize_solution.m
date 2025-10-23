function fig = visualize_solution(u, xx, yy, nblock, nnodes, options)
%VISUALIZE_SOLUTION 可视化位移场和变形（使用 contourf 绘制）
%
% fig = visualize_solution(u, xx, yy, nblock, nnodes, options)
%
% INPUTS
%   u       : (2*nblock*nnodes) × 1 位移向量 [u1,v1,u2,v2,...,un,vn]
%   xx      : cell 数组，每个单元的配点 x 坐标（nt×nt 矩阵）
%   yy      : cell 数组，每个单元的配点 y 坐标（nt×nt 矩阵）
%   nblock  : 单元数量
%   nnodes  : 每个单元的配点数
%   options : 可选参数结构体
%             .scale : 变形放大系数（默认自动）
%             .showUndeformed : 是否显示未变形构型（默认 true）
%             .nlevels : 等高线层数（默认 20）
%
% OUTPUTS
%   fig : 单一 figure 句柄（包含 2×2 子图）
%
% DESCRIPTION
%   使用 contourf 绘制云图，所有块体共享同一个 colorbar。
%   子图布局：[u方向位移, v方向位移; 位移幅值, 变形构型]
%
% EXAMPLE
%   fig = visualize_solution(u, xx, yy, nblock, nnodes);
%
% Author: (W. Huang)
% -------------------------------------------------------------------------

if nargin < 6
    options = struct();
end

% 默认选项
if ~isfield(options, 'showUndeformed'), options.showUndeformed = true; end
if ~isfield(options, 'nlevels'), options.nlevels = 20; end

% 提取 u 和 v 分量
u_disp = u(1:2:end);  % u 方向位移
v_disp = u(2:2:end);  % v 方向位移

% 计算每个方向的配点数
nt = sqrt(nnodes);

% 确保 nt 是整数
if abs(nt - round(nt)) > 1e-10
    error('nnodes (%d) 不是完全平方数！nt = %.4f', nnodes, nt);
end
nt = round(nt);

% 重排坐标和位移数据为 nt×nt 矩阵
% 注意：xx{i} 和 yy{i} 从 mappings 返回的是列向量，需要 reshape 成矩阵
xx_mat = cell(1, nblock);
yy_mat = cell(1, nblock);
u_blocks = cell(1, nblock);
v_blocks = cell(1, nblock);

for i = 1:nblock
    % 坐标数据 reshape 成 nt×nt 矩阵（用于 contourf）
    xx_mat{i} = reshape(xx{i}, [nt, nt]);
    yy_mat{i} = reshape(yy{i}, [nt, nt]);
    
    % 位移数据也 reshape 成 nt×nt 矩阵
    idx_start = (i-1)*nnodes + 1;
    idx_end = i*nnodes;
    u_blocks{i} = reshape(u_disp(idx_start:idx_end), [nt, nt]);
    v_blocks{i} = reshape(v_disp(idx_start:idx_end), [nt, nt]);
end

% 计算位移幅值
disp_magnitude = sqrt(u_disp.^2 + v_disp.^2);
disp_mag_blocks = cell(1, nblock);
for i = 1:nblock
    idx_start = (i-1)*nnodes + 1;
    idx_end = i*nnodes;
    disp_mag_blocks{i} = reshape(disp_magnitude(idx_start:idx_end), [nt, nt]);
end

% 自动确定变形放大系数
if ~isfield(options, 'scale')
    max_disp = max(abs(u));
    xx_all = cell2mat(cellfun(@(x) x(:), xx_mat, 'UniformOutput', false));
    yy_all = cell2mat(cellfun(@(y) y(:), yy_mat, 'UniformOutput', false));
    domain_size = max([range(xx_all), range(yy_all)]);
    if max_disp > 0
        options.scale = domain_size / max_disp * 0.1;  % 10% 的域尺寸
    else
        options.scale = 1.0;
    end
end

% 计算全局颜色范围（确保所有子图使用相同的颜色映射）
u_min = min(u_disp);
u_max = max(u_disp);
v_min = min(v_disp);
v_max = max(v_disp);
mag_min = min(disp_magnitude);
mag_max = max(disp_magnitude);

%% 创建单一图形窗口，包含 2×2 子图
fig = figure('Name', '位移场可视化', 'Position', [100, 100, 1400, 1000]/2);

%% 子图 1: u 方向位移云图
subplot(2, 2, 1);
hold on;
for i = 1:nblock
    % 使用 contourf 绘制填充等高线
    contourf(xx_mat{i}, yy_mat{i}, u_blocks{i}, options.nlevels, 'LineStyle', 'none');
end
colormap('jet');
caxis([u_min, u_max]);
cb1 = colorbar;
cb1.Label.String = 'u 位移';
cb1.Label.FontSize = 11;
title(sprintf('u 方向位移 (最大: %.3e)', max(abs(u_disp))), 'FontSize', 13, 'FontWeight', 'bold');
xlabel('X', 'FontSize', 11);
ylabel('Y', 'FontSize', 11);
axis equal tight;
grid on;
hold off;

%% 子图 2: v 方向位移云图
subplot(2, 2, 2);
hold on;
for i = 1:nblock
    contourf(xx_mat{i}, yy_mat{i}, v_blocks{i}, options.nlevels, 'LineStyle', 'none');
end
colormap('jet');
caxis([v_min, v_max]);
cb2 = colorbar;
cb2.Label.String = 'v 位移';
cb2.Label.FontSize = 11;
title(sprintf('v 方向位移 (最大: %.3e)', max(abs(v_disp))), 'FontSize', 13, 'FontWeight', 'bold');
xlabel('X', 'FontSize', 11);
ylabel('Y', 'FontSize', 11);
axis equal tight;
grid on;
hold off;

%% 子图 3: 位移幅值云图
subplot(2, 2, 3);
hold on;
for i = 1:nblock
    contourf(xx_mat{i}, yy_mat{i}, disp_mag_blocks{i}, options.nlevels, 'LineStyle', 'none');
end
colormap('jet');
caxis([mag_min, mag_max]);
cb3 = colorbar;
cb3.Label.String = '|u|';
cb3.Label.FontSize = 11;
title(sprintf('位移幅值 (最大: %.3e)', max(disp_magnitude)), 'FontSize', 13, 'FontWeight', 'bold');
xlabel('X', 'FontSize', 11);
ylabel('Y', 'FontSize', 11);
axis equal tight;
grid on;
hold off;

%% 子图 4: 变形构型
subplot(2, 2, 4);
hold on;

% 未变形构型（使用浅灰色线框）
if options.showUndeformed
    for i = 1:nblock
        plot(xx_mat{i}, yy_mat{i}, 'Color', [0.7, 0.7, 0.7], 'LineWidth', 0.5);
        plot(xx_mat{i}', yy_mat{i}', 'Color', [0.7, 0.7, 0.7], 'LineWidth', 0.5);
    end
end

% 变形后构型（使用 contourf 显示位移幅值）
for i = 1:nblock
    % 计算变形后的坐标
    xx_def = xx_mat{i} + u_blocks{i} * options.scale;
    yy_def = yy_mat{i} + v_blocks{i} * options.scale;
    
    % 绘制变形后的云图
    contourf(xx_def, yy_def, disp_mag_blocks{i}, options.nlevels, 'LineStyle', 'none');
end

colormap('jet');
caxis([mag_min, mag_max]);
cb4 = colorbar;
cb4.Label.String = '|u|';
cb4.Label.FontSize = 11;
title(sprintf('变形构型 (放大: %.1f×)', options.scale), 'FontSize', 13, 'FontWeight', 'bold');
xlabel('X', 'FontSize', 11);
ylabel('Y', 'FontSize', 11);
axis equal tight;
grid on;
hold off;

% 添加总标题
sgtitle('位移场和变形可视化', 'FontSize', 16, 'FontWeight', 'bold');

fprintf('  生成位移场可视化图形（4个子图）\n');

end
