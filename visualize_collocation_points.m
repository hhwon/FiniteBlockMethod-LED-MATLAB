function fig = visualize_collocation_points(xx, yy, options)
%VISUALIZE_COLLOCATION_POINTS 可视化各单元的配点分布
%
% fig = visualize_collocation_points(xx, yy, options)
%
% INPUTS
%   xx      : cell 数组，每个单元的配点 x 坐标
%   yy      : cell 数组，每个单元的配点 y 坐标
%   options : 可选参数结构体
%             .markerSize : 配点标记大小（默认 4）
%             .showGrid   : 是否显示网格线（默认 false）
%             .figTitle   : 图形标题（默认自动生成）
%             .figName    : 图形窗口名称（默认 '配点分布'）
%
% OUTPUTS
%   fig : figure 句柄
%
% EXAMPLE
%   [xx, yy, ~, ~, ~] = generate_collocation_points(xcoor, ycoor, type, nt);
%   visualize_collocation_points(xx, yy);
%
% See also: generate_collocation_points
%
% Author: (W. Huang)
% -------------------------------------------------------------------------

%% 参数处理
if nargin < 3
    options = struct();
end

% 默认选项
if ~isfield(options, 'markerSize'), options.markerSize = 4; end
if ~isfield(options, 'showGrid'),   options.showGrid   = false; end
if ~isfield(options, 'figName'),    options.figName    = '配点分布'; end

nblock = numel(xx);

if ~isfield(options, 'figTitle')
    total_points = sum(cellfun(@numel, xx));
    options.figTitle = sprintf('配点分布 (共 %d 个块, %d 个配点)', nblock, total_points);
end

%% 创建图形
fig = figure('Name', options.figName);
hold on;

% 生成颜色映射
colors = lines(nblock);

% 绘制每个块的配点
legend_entries = cell(1, nblock);
for i = 1:nblock
    npts = numel(xx{i});
    
    if options.showGrid
        % 显示网格线
        plot(xx{i}, yy{i}, '.', ...
             'Color', colors(i,:), ...
             'MarkerSize', options.markerSize, ...
             'LineWidth', 0.5);
    else
        % 只显示配点
        plot(xx{i}(:), yy{i}(:), '.', ...
             'Color', colors(i,:), ...
             'MarkerSize', options.markerSize);
    end
    
    legend_entries{i} = sprintf('Block %d (%d points)', i, npts);
end

%% 图形美化
axis equal;
grid on;
xlabel('X', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y', 'FontSize', 12, 'FontWeight', 'bold');
title(options.figTitle, 'FontSize', 14, 'FontWeight', 'bold');
legend(legend_entries, 'Location', 'best', 'FontSize', 10);
hold off;

end
