function [Aphi, ADphix, ADphiy, ADphixx, ADphixy, ADphiyy] = ...
    assemble_global_operators(Phi, Dphix, Dphiy, Dphixx, Dphixy, Dphiyy)
%ASSEMBLE_GLOBAL_OPERATORS 装配全局形函数算子矩阵
%
% [Aphi, ADphix, ADphiy, ADphixx, ADphixy, ADphiyy] = ...
%     assemble_global_operators(Phi, Dphix, Dphiy, Dphixx, Dphixy, Dphiyy)
%
% INPUTS
%   Phi    : cell 数组，每个单元的基函数值矩阵
%   Dphix  : cell 数组，每个单元的 x 方向一阶导数算子
%   Dphiy  : cell 数组，每个单元的 y 方向一阶导数算子
%   Dphixx : cell 数组，每个单元的 x 方向二阶导数算子
%   Dphixy : cell 数组，每个单元的混合二阶导数算子
%   Dphiyy : cell 数组，每个单元的 y 方向二阶导数算子
%
% OUTPUTS
%   Aphi   : 全局基函数值矩阵（块对角）
%   ADphix : 全局 x 方向一阶导数算子（块对角）
%   ADphiy : 全局 y 方向一阶导数算子（块对角）
%   ADphixx: 全局 x 方向二阶导数算子（块对角）
%   ADphixy: 全局混合二阶导数算子（块对角）
%   ADphiyy: 全局 y 方向二阶导数算子（块对角）
%
% DESCRIPTION
%   本函数使用块对角矩阵（blkdiag）将各单元的局部算子装配为全局算子。
%   每个单元的自由度独立编号，全局矩阵维度为 (nblock*nnodes) × (nblock*nnodes)。
%
%   装配后的全局算子可用于：
%   - 刚度矩阵装配
%   - 边界条件施加
%   - 求解线性方程组
%
% EXAMPLE
%   [Phi, Dx, Dy, Dxx, Dxy, Dyy, ~, ~, ~, ~] = ...
%       compute_shape_functions(xcoor, ycoor, type, xt, yt, span);
%   [Aphi, ADx, ADy, ADxx, ADxy, ADyy] = ...
%       assemble_global_operators(Phi, Dx, Dy, Dxx, Dxy, Dyy);
%
% See also: blkdiag, compute_shape_functions
%
% Author: (W. Huang)
% -------------------------------------------------------------------------

%% 使用块对角矩阵装配全局算子
Aphi    = blkdiag(Phi{:});
ADphix  = blkdiag(Dphix{:});
ADphiy  = blkdiag(Dphiy{:});
ADphixx = blkdiag(Dphixx{:});
ADphixy = blkdiag(Dphixy{:});
ADphiyy = blkdiag(Dphiyy{:});

end
