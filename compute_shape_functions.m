function [Phi, Dphix, Dphiy, Dphixx, Dphixy, Dphiyy, delta11, delta12, delta21, delta22] = ...
    compute_shape_functions(xcoor, ycoor, type, xt, yt, span)
%COMPUTE_SHAPE_FUNCTIONS 计算所有单元的形函数及其导数
%
% [Phi, Dphix, Dphiy, Dphixx, Dphixy, Dphiyy, delta11, delta12, delta21, delta22] = ...
%     compute_shape_functions(xcoor, ycoor, type, xt, yt, span)
%
% INPUTS
%   xcoor : cell 数组，每个单元的节点 x 坐标
%   ycoor : cell 数组，每个单元的节点 y 坐标
%   type  : cell 数组，单元类型标识
%   xt    : nt×nt 矩阵，参考空间 ξ 坐标网格
%   yt    : nt×nt 矩阵，参考空间 η 坐标网格
%   span  : 行向量，一维离散化节点
%
% OUTPUTS
%   Phi    : cell 数组，每个单元的基函数值矩阵 (nnodes×nnodes)
%   Dphix  : cell 数组，x 方向一阶导数算子 (nnodes×nnodes)
%   Dphiy  : cell 数组，y 方向一阶导数算子 (nnodes×nnodes)
%   Dphixx : cell 数组，x 方向二阶导数算子 (nnodes×nnodes)
%   Dphixy : cell 数组，混合二阶导数算子 (nnodes×nnodes)
%   Dphiyy : cell 数组，y 方向二阶导数算子 (nnodes×nnodes)
%   delta11: cell 数组，逆雅可比分量 ∂ξ/∂x (nnodes×nnodes)
%   delta12: cell 数组，逆雅可比分量 ∂η/∂x (nnodes×nnodes)
%   delta21: cell 数组，逆雅可比分量 ∂ξ/∂y (nnodes×nnodes)
%   delta22: cell 数组，逆雅可比分量 ∂η/∂y (nnodes×nnodes)
%
% DESCRIPTION
%   本函数执行以下步骤：
%   1. 对每个单元，计算等参映射的雅可比矩阵（通过 mappings）
%   2. 使用雅可比信息，计算形函数及其导数（通过 diffLarg）
%   3. 返回所有单元的形函数算子，供后续刚度矩阵装配使用
%
%   形函数计算基于张量积 Lagrange 基函数：
%   - φ(x,y) = Lx(ξ) ⊗ Ly(η)
%   - 通过链式法则转换到物理坐标系
%
% EXAMPLE
%   [xcoor, ycoor, type] = expand_nodes(CASE.blocks);
%   [xx, yy, xt, yt, span] = generate_collocation_points(xcoor, ycoor, type, CASE.nt);
%   [Phi, Dx, Dy, Dxx, Dxy, Dyy, d11, d12, d21, d22] = ...
%       compute_shape_functions(xcoor, ycoor, type, xt, yt, span);
%
% See also: mappings, diffLarg, larg
%
% Author: (W. Huang)
% -------------------------------------------------------------------------

nblock = numel(xcoor);

% 预分配输出 cell 数组
delta11 = cell(1, nblock);
delta12 = cell(1, nblock);
delta21 = cell(1, nblock);
delta22 = cell(1, nblock);
Phi     = cell(1, nblock);
Dphix   = cell(1, nblock);
Dphiy   = cell(1, nblock);
Dphixx  = cell(1, nblock);
Dphixy  = cell(1, nblock);
Dphiyy  = cell(1, nblock);

%% 步骤 1: 计算每个单元的雅可比矩阵（等参映射）
fprintf('  步骤 1/2: 计算雅可比矩阵...\n');
for i = 1:nblock
    % 调用 mappings 函数获取逆雅可比分量
    [delta11{i}, delta12{i}, delta21{i}, delta22{i}, ~, ~] = ...
        mappings(xcoor{i}, ycoor{i}, xt, yt, type{i});
end

%% 步骤 2: 计算形函数及其导数
fprintf('  步骤 2/2: 计算形函数及导数...\n');
for i = 1:nblock
    % 调用 diffLarg 计算基函数及其导数算子
    [Phi{i}, Dphix{i}, Dphiy{i}, Dphixx{i}, Dphixy{i}, Dphiyy{i}] = ...
        diffLarg(xt, yt, span, delta11{i}, delta12{i}, delta21{i}, delta22{i});
end

end
