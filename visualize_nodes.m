function fig = visualize_nodes(xcoor, ycoor, type, options)
%VISUALIZE_NODES 可视化扩展后的节点和单元结构
%
% fig = visualize_nodes(xcoor, ycoor, type, options)
%
% INPUTS
%   xcoor   : cell 数组，每个单元的 x 坐标
%   ycoor   : cell 数组，每个单元的 y 坐标
%   type    : cell 数组，单元类型标识
%   options : 可选参数结构体
%             .showNodeLabels  : 是否显示节点编号 (默认 true)
%             .showBlockLabels : 是否显示块标签 (默认 true)
%             .figTitle        : 图形标题 (默认自动生成)
%             .figName         : 图形窗口名称 (默认 '节点可视化')
%
% OUTPUTS
%   fig : figure 句柄
%
% EXAMPLE
%   [xcoor, ycoor, type] = expand_nodes(CASE.blocks);
%   visualize_nodes(xcoor, ycoor, type);
%
% See also: expand_nodes, plot_node_summary
%
% Author: (W. Huang)
% -------------------------------------------------------------------------

%% 参数处理
if nargin < 4
    options = struct();
end

% 默认选项
if ~isfield(options, 'showNodeLabels'),  options.showNodeLabels  = true; end
if ~isfield(options, 'showBlockLabels'), options.showBlockLabels = true; end
if ~isfield(options, 'figName'),         options.figName         = '节点可视化'; end

nblock = numel(xcoor);

if ~isfield(options, 'figTitle')
    options.figTitle = sprintf('节点扩展结果 (共 %d 个块)', nblock);
end

%% 创建图形
fig = figure('Name', options.figName);
hold on;

% 生成颜色映射
colors = lines(nblock);

% 绘制每个块
legend_entries = cell(1, nblock);
for i = 1:nblock
    % 绘制节点和连线
    plot(xcoor{i}, ycoor{i}, 'o', ...
         'Color', colors(i,:), ...
         'MarkerSize', 8, ...
         'LineWidth', 1.5, ...
         'MarkerFaceColor', colors(i,:), ...
         'DisplayName', sprintf('Block %d (%d nodes)', i, numel(xcoor{i})));
    
    % 标注节点编号
    if options.showNodeLabels
        for j = 1:numel(xcoor{i})
            text(xcoor{i}(j), ycoor{i}(j), sprintf(' %d', j), ...
                 'FontSize', 15, ...
                 'Color', colors(i,:), ...
                 'FontWeight', 'bold');
        end
    end
    
    % 标注块信息
    if options.showBlockLabels
        xc = mean(xcoor{i});
        yc = mean(ycoor{i});
        text(xc, yc, sprintf('Block %d\n(type=%d)', i, type{i}), ...
             'FontSize', 11, ...
             'HorizontalAlignment', 'center', ...
             'BackgroundColor', 'white', ...
             'EdgeColor', colors(i,:), ...
             'Margin', 2);
    end
    
    legend_entries{i} = sprintf('Block %d (%d nodes)', i, numel(xcoor{i}));
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
