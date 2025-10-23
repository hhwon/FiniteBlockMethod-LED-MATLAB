function fig = visualize_operator_sparsity(Aphi, ADphix, ADphiy, options)
%VISUALIZE_OPERATOR_SPARSITY 可视化全局算子矩阵的稀疏结构
%
% fig = visualize_operator_sparsity(Aphi, ADphix, ADphiy, options)
%
% INPUTS
%   Aphi    : 全局基函数值矩阵
%   ADphix  : 全局 x 方向一阶导数算子
%   ADphiy  : 全局 y 方向一阶导数算子
%   options : 可选参数结构体
%             .showColorbar : 是否显示色条（默认 true）
%             .figName      : 图形窗口名称（默认 '算子稀疏结构'）
%
% OUTPUTS
%   fig : figure 句柄
%
% EXAMPLE
%   [Aphi, ADx, ADy, ~, ~, ~] = assemble_global_operators(...);
%   visualize_operator_sparsity(Aphi, ADx, ADy);
%
% See also: spy, assemble_global_operators
%
% Author: (W. Huang)
% -------------------------------------------------------------------------

%% 参数处理
if nargin < 4
    options = struct();
end

if ~isfield(options, 'showColorbar'), options.showColorbar = true; end
if ~isfield(options, 'figName'),      options.figName      = '算子稀疏结构'; end

%% 创建图形
fig = figure('Name', options.figName, 'Position', [100, 100, 1200, 400]);

%% 子图 1: 基函数值矩阵 Aphi
subplot(1, 3, 1);
spy(Aphi);
title(sprintf('Φ 基函数值矩阵\n非零元素: %d (%.2f%%)', ...
    nnz(Aphi), 100*nnz(Aphi)/numel(Aphi)), 'FontSize', 12);
xlabel('列索引', 'FontSize', 10);
ylabel('行索引', 'FontSize', 10);
grid on;

%% 子图 2: x 方向导数算子 ADphix
subplot(1, 3, 2);
spy(ADphix);
title(sprintf('∂Φ/∂x 导数算子\n非零元素: %d (%.2f%%)', ...
    nnz(ADphix), 100*nnz(ADphix)/numel(ADphix)), 'FontSize', 12);
xlabel('列索引', 'FontSize', 10);
ylabel('行索引', 'FontSize', 10);
grid on;

%% 子图 3: y 方向导数算子 ADphiy
subplot(1, 3, 3);
spy(ADphiy);
title(sprintf('∂Φ/∂y 导数算子\n非零元素: %d (%.2f%%)', ...
    nnz(ADphiy), 100*nnz(ADphiy)/numel(ADphiy)), 'FontSize', 12);
xlabel('列索引', 'FontSize', 10);
ylabel('行索引', 'FontSize', 10);
grid on;

%% 总标题
sgtitle('全局算子稀疏结构分析', 'FontSize', 14, 'FontWeight', 'bold');

end
