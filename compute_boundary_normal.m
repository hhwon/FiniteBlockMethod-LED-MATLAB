function [nx, ny] = compute_boundary_normal(xi, eta, delta11, delta12, delta21, delta22, xt, yt, face)
%COMPUTE_BOUNDARY_NORMAL_V2 计算边界上任意点的外法向量（使用已有的逆雅可比矩阵）
%
% [nx, ny] = compute_boundary_normal_v2(xi, eta, delta11, delta12, delta21, delta22, xt, yt, face)
%
% INPUTS
%   xi, eta : 标量，参数空间坐标 (∈ [-1, 1])
%   delta11 : 矩阵 (nt×nt)，逆雅可比分量 ∂ξ/∂x
%   delta12 : 矩阵 (nt×nt)，逆雅可比分量 ∂η/∂x
%   delta21 : 矩阵 (nt×nt)，逆雅可比分量 ∂ξ/∂y
%   delta22 : 矩阵 (nt×nt)，逆雅可比分量 ∂η/∂y
%   xt, yt  : 矩阵 (nt×nt)，参数空间的网格坐标
%   face    : 字符，边界面标识 ('L', 'R', 'B', 'T')
%             'L' = 左边界 (xi = -1)
%             'R' = 右边界 (xi = 1)
%             'B' = 下边界 (eta = -1)
%             'T' = 上边界 (eta = 1)
%
% OUTPUTS
%   nx, ny : 标量，归一化的外法向量分量
%
% DESCRIPTION
%   本函数直接使用 compute_shape_functions 计算的逆雅可比矩阵，避免重复插值。
%   
%   核心思想：
%   1. 从逆雅可比矩阵反推正雅可比矩阵：
%      逆雅可比矩阵 J^{-1} = [delta11  delta12]
%                              [delta21  delta22]
%      
%      正雅可比矩阵 J = [∂x/∂ξ  ∂y/∂ξ] = (1/det) * [ delta22  -delta12]
%                       [∂x/∂η  ∂y/∂η]              [-delta21   delta11]
%      
%      其中 det = delta11*delta22 - delta12*delta21
%   
%   2. 提取切向量：
%      - 左/右边界 (xi = ±1)：切向 t = [∂x/∂η, ∂y/∂η]
%      - 下/上边界 (eta = ±1)：切向 t = [∂x/∂ξ, ∂y/∂ξ]
%   
%   3. 旋转 90°得到法向量并归一化：
%      - 左边界：n = [-ty, tx] (指向左侧 -x)
%      - 右边界：n = [ty, -tx] (指向右侧 +x)
%      - 下边界：n = [ty, -tx] (指向下方 -y)
%      - 上边界：n = [-ty, tx] (指向上方 +y)
%
% NOTES
%   - 使用最近邻插值在 delta 矩阵中查找对应 (xi, eta) 的值
%   - 假设 xt, yt 矩阵对应的参数网格在 [-1, 1]² 内
%   - 相比 compute_boundary_normal，本版本避免重复计算，更高效且一致
%
% EXAMPLE
%   % 假设已经调用过 compute_shape_functions
%   [~, ~, ~, ~, ~, ~, delta11, delta12, delta21, delta22] = ...
%       compute_shape_functions(xcoor, ycoor, type, xt, yt, span);
%   
%   % 计算左边界中点的法向量
%   xi = -1; eta = 0;
%   [nx, ny] = compute_boundary_normal_v2(xi, eta, delta11{iblock}, delta12{iblock}, ...
%                                         delta21{iblock}, delta22{iblock}, xt, yt, 'L');
%
% See also: compute_shape_functions, mappings, assemble_global_stiffness
%
% Author: (W. Huang)
% Date: 2025-10-17
% -------------------------------------------------------------------------

%% 1. 在网格中找到最接近 (xi, eta) 的点

% 计算与所有网格点的距离
dist = (xt - xi).^2 + (yt - eta).^2;

% 找到最近的点的索引
[~, idx] = min(dist(:));

%% 2. 提取该点的逆雅可比矩阵分量

d11 = delta11(idx);  % ∂ξ/∂x
d12 = delta12(idx);  % ∂η/∂x
d21 = delta21(idx);  % ∂ξ/∂y
d22 = delta22(idx);  % ∂η/∂y

%% 3. 反推正雅可比矩阵

% 计算行列式
det_J_inv = d11 * d22 - d12 * d21;

% 正雅可比矩阵分量
% J = [∂x/∂ξ  ∂y/∂ξ]
%     [∂x/∂η  ∂y/∂η]
dxdxi  =  d22 / det_J_inv;   % ∂x/∂ξ
dydxi  = -d12 / det_J_inv;   % ∂y/∂ξ
dxdeta = -d21 / det_J_inv;   % ∂x/∂η
dydeta =  d11 / det_J_inv;   % ∂y/∂η

%% 4. 根据边界位置提取切向量并计算法向量

switch upper(face)
    case 'L'  % 左边界 (xi = -1)
        % 切向量沿 η 方向：t = [∂x/∂η, ∂y/∂η]
        tx = dxdeta;
        ty = dydeta;
        % 法向量：旋转90°并指向左侧（-x方向）
        nx = -ty;
        ny = tx;
        
    case 'R'  % 右边界 (xi = 1)
        % 切向量沿 η 方向：t = [∂x/∂η, ∂y/∂η]
        tx = dxdeta;
        ty = dydeta;
        % 法向量：旋转-90°并指向右侧（+x方向）
        nx = ty;
        ny = -tx;
        
    case 'B'  % 下边界 (eta = -1)
        % 切向量沿 ξ 方向：t = [∂x/∂ξ, ∂y/∂ξ]
        tx = dxdxi;
        ty = dydxi;
        % 法向量：旋转-90°并指向下方（-y方向）
        nx = ty;
        ny = -tx;
        
    case 'T'  % 上边界 (eta = 1)
        % 切向量沿 ξ 方向：t = [∂x/∂ξ, ∂y/∂ξ]
        tx = dxdxi;
        ty = dydxi;
        % 法向量：旋转90°并指向上方（+y方向）
        nx = -ty;
        ny = tx;
        
    otherwise
        error('compute_boundary_normal_v2: 未知的边界标识 face = %s （应为 L/R/B/T）', face);
end

%% 5. 归一化

norm_n = sqrt(nx^2 + ny^2);
if norm_n < 1e-12
    warning('compute_boundary_normal_v2: 法向量长度接近零 (%.2e)，可能存在退化映射', norm_n);
    norm_n = 1;  % 避免除以零
end

nx = nx / norm_n;
ny = ny / norm_n;

end
