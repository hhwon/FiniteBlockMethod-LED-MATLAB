function [xx, yy, xt, yt, span] = generate_collocation_points(xcoor, ycoor, type, nt)
%GENERATE_COLLOCATION_POINTS 在各单元内生成配点（collocation points）
%
% [xx, yy, xt, yt, span] = generate_collocation_points(xcoor, ycoor, type, nt)
%
% INPUTS
%   xcoor : cell 数组，每个单元的节点 x 坐标（已扩展）
%   ycoor : cell 数组，每个单元的节点 y 坐标（已扩展）
%   type  : cell 数组，单元类型标识
%   nt    : 标量，每个方向的离散化节点数（从 CASE.nt 获取）
%
% OUTPUTS
%   xx    : cell 数组，每个单元的物理空间配点 x 坐标（nt×nt 矩阵）
%   yy    : cell 数组，每个单元的物理空间配点 y 坐标（nt×nt 矩阵）
%   xt    : nt×nt 矩阵，参考空间 ξ 坐标网格（所有单元共享）
%   yt    : nt×nt 矩阵，参考空间 η 坐标网格（所有单元共享）
%   span  : 行向量，[-1,1] 区间的离散化节点（采用正弦分布）
%
% DESCRIPTION
%   本函数通过以下步骤生成配点：
%   1. 在参考空间 [-1,1]×[-1,1] 上生成 nt×nt 网格（正弦分布）
%   2. 对每个单元，通过等参映射 mappings() 将参考坐标映射到物理空间
%   3. 返回物理空间坐标和参考空间坐标供后续计算使用
%
% EXAMPLE
%   CASE = load_json_case('case_linear.json');
%   [xcoor, ycoor, type] = expand_nodes(CASE.blocks);
%   [xx, yy, xt, yt, span] = generate_collocation_points(xcoor, ycoor, type, CASE.nt);
%
% See also: mappings, expand_nodes
%
% Author: (W. Huang)
% -------------------------------------------------------------------------

%% 生成参考空间网格（正弦分布，端点密集）
theta = linspace(-pi/2, pi/2, nt);
span = sin(theta);  % 行向量，长度 nt，范围 [-1, 1]

% 生成 2D 网格
[xt, yt] = meshgrid(span, span);  % 两者均为 nt×nt 矩阵

%% 对每个单元执行等参映射
nblock = numel(xcoor);
xx = cell(1, nblock);
yy = cell(1, nblock);

for i = 1:nblock
    % 调用 mappings 函数进行等参映射
    % mappings 返回：[delta11, delta12, delta21, delta22, xx_mapped, yy_mapped]
    [~, ~, ~, ~, xx{i}, yy{i}] = mappings(xcoor{i}, ycoor{i}, xt, yt, type{i});
end

end
