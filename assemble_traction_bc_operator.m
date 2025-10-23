function k = assemble_traction_bc_operator(phi, dphix, dphiy, x, ndof, normal, E, nu, plane_strain, axisymmetric, r)
%ASSEMBLE_TRACTION_BC_OPERATOR 装配边界牵引力条件的算子行
%
% k = assemble_traction_bc_operator(phi, dphix, dphiy, x, ndof, normal, E, nu, plane_strain, axisymmetric, r)
%
% INPUTS
%   phi    : 1×n_nodes 行向量，配点处的基函数值
%   dphix  : 1×n_nodes 行向量，配点处 ∂φ/∂x (轴对称时为 ∂φ/∂r)
%   dphiy  : 1×n_nodes 行向量，配点处 ∂φ/∂y (轴对称时为 ∂φ/∂z)
%   x      : 标量，配点的 x 坐标（可用于变系数，此处未用）
%   ndof   : 标量，节点总数（注意：全局自由度数 = 2*ndof）
%   normal : 2×1 向量，边界外法向量 [nx; ny] 或 [nr; nz]
%   E      : 标量，杨氏模量
%   nu     : 标量，泊松比
%   plane_strain : 逻辑值，是否为平面应变（默认 true）
%   axisymmetric : 逻辑值，是否为轴对称问题（默认 false）
%   r      : 标量，配点的径向坐标（仅轴对称时使用，必须 > 0）
%
% OUTPUTS
%   k : 2×(2*ndof) 矩阵，对应该配点的牵引力边界条件
%       列按交错顺序排列：[u1, v1, u2, v2, ..., un, vn] 或 [ur1, uz1, ur2, uz2, ...]
%
% DESCRIPTION
%   本函数实现边界牵引力条件：σ·n = t_prescribed
%
%   平面应变/应力：
%     应力：σxx = (λ+2μ)εxx + λεyy, σyy = (λ+2μ)εyy + λεxx, σxy = 2μεxy
%     应变：εxx = ∂u/∂x, εyy = ∂v/∂y, εxy = (∂u/∂y + ∂v/∂x)/2
%     牵引力：tx = σxx*nx + σxy*ny, ty = σxy*nx + σyy*ny
%
%   轴对称（r-z 坐标）：
%     应力：σrr = (λ+2μ)εrr + λεzz + λεθθ, σzz = λεrr + (λ+2μ)εzz + λεθθ
%           σrz = 2μεrz, σθθ = λεrr + λεzz + (λ+2μ)εθθ
%     应变：εrr = ∂u/∂r, εzz = ∂w/∂z, εθθ = u/r, εrz = (∂u/∂z + ∂w/∂r)/2
%     牵引力：tr = σrr*nr + σrz*nz, tz = σrz*nr + σzz*nz
%
%   注意：当 r→0 时，使用 L'Hospital 法则：u/r → ∂u/∂r (即 dphix)
%
% See also: assemble_strong_form_operator
%
% Author: (W. Huang)
% -------------------------------------------------------------------------

if nargin < 9
    plane_strain = true;
end
if nargin < 10
    axisymmetric = false;
end
if nargin < 11
    r = 1.0;  % 默认值（平面问题时不使用）
end

% 计算 Lamé 常数
if plane_strain
    lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
    mu = E / (2 * (1 + nu));
else
    % 平面应力
    lambda = E * nu / (1 - nu^2);
    mu = E / (2 * (1 + nu));
end

% 提取法向分量
% normal 是 2×1 向量 [nx; ny] 或 [nr; nz]
% 在平面问题中：nx, ny 表示 x, y 方向的法向量分量
% 在轴对称问题中：nx→nr, ny→nz 表示 r, z 方向的法向量分量
nx = normal(1);  % x 方向（平面）或 r 方向（轴对称）
ny = normal(2);  % y 方向（平面）或 z 方向（轴对称）

% 注意：ndof 是节点数，但全局自由度是 2*ndof（交错排列 u1,v1,u2,v2,...）
n_nodes = ndof;
n_dof_global = 2 * n_nodes;

% 预分配刚度矩阵行（2 行，每行对应一个方向，列数为全局自由度数）
k = zeros(2, n_dof_global);

if axisymmetric
    % ========================================================================
    % 轴对称边界牵引力（r-z 坐标系，u=径向位移，w=轴向位移）
    % ========================================================================
    % 应力分量：
    %   σrr = (λ+2μ)∂u/∂r + λ∂w/∂z + λu/r
    %   σzz = λ∂u/∂r + (λ+2μ)∂w/∂z + λu/r
    %   σrz = μ(∂u/∂z + ∂w/∂r)
    %   σθθ = λ∂u/∂r + λ∂w/∂z + (λ+2μ)u/r
    %
    % 牵引力：
    %   tr = σrr*nr + σrz*nz
    %   tz = σrz*nr + σzz*nz
    %
    % 特殊处理：当 r→0 时，u/r → ∂u/∂r (即 dphix)
    % ========================================================================
    
    % 在轴对称情况下重命名以提高可读性
    nr = nx;  % 径向法向量分量
    nz = ny;  % 轴向法向量分量
    
    % 数值稳定性阈值
    r_tol = 1e-5;
    
    if abs(r) < r_tol
        % 对称轴上（r ≈ 0）的极限情况
        % εθθ = u/r → ∂u/∂r
        % σrr = (λ+2μ)∂u/∂r + λ∂w/∂z + λ∂u/∂r = (λ+3μ+λ)∂u/∂r + λ∂w/∂z
        % σzz = λ∂u/∂r + (λ+2μ)∂w/∂z + λ∂u/∂r = 2λ∂u/∂r + (λ+2μ)∂w/∂z
        
        % r 方向牵引力: tr = σrr*nr + σrz*nz
        % 对 u 的贡献（奇数列）
        k(1, 1:2:end) = ((lambda+2*mu)*dphix + lambda*dphix)*nr + mu*dphiy*nz;
        % 对 w 的贡献（偶数列）
        k(1, 2:2:end) = lambda*dphiy*nr + mu*dphix*nz;
        
        % z 方向牵引力: tz = σrz*nr + σzz*nz
        % 对 u 的贡献（奇数列）
        k(2, 1:2:end) = mu*dphiy*nr + (2*lambda*dphix)*nz;
        % 对 w 的贡献（偶数列）
        k(2, 2:2:end) = mu*dphix*nr + (lambda+2*mu)*dphiy*nz;
    else
        % 一般情况（r > 0）
        % r 方向牵引力: tr = σrr*nr + σrz*nz
        % 对 u 的贡献（奇数列）
        k(1, 1:2:end) = ((lambda+2*mu)*dphix + lambda*phi/r)*nr + mu*dphiy*nz;
        % 对 w 的贡献（偶数列）
        k(1, 2:2:end) = lambda*dphiy*nr + mu*dphix*nz;
        
        % z 方向牵引力: tz = σrz*nr + σzz*nz
        % 对 u 的贡献（奇数列）
        k(2, 1:2:end) = mu*dphiy*nr + (lambda*dphix + lambda*phi/r)*nz;
        % 对 w 的贡献（偶数列）
        k(2, 2:2:end) = mu*dphix*nr + (lambda+2*mu)*dphiy*nz;
        % 对 w 的贡献（偶数列）
        k(2, 2:2:end) = mu*dphix*nx + (lambda+2*mu)*dphiy*nz;
    end
else
    % ========================================================================
    % 平面应变/平面应力（x-y 坐标系，u=x方向位移，v=y方向位移）
    % ========================================================================
    % x 方向牵引力: tx = σxx*nx + σxy*ny
    % σxx = (λ+2μ)∂u/∂x + λ∂v/∂y
    % σxy = μ(∂u/∂y + ∂v/∂x)
    
    % 对 u 的贡献（奇数列：1,3,5,...）
    k(1, 1:2:end) = ((lambda + 2*mu) * dphix * nx + mu * dphiy * ny);
    % 对 v 的贡献（偶数列：2,4,6,...）
    k(1, 2:2:end) = (lambda * dphiy * nx + mu * dphix * ny);
    
    % y 方向牵引力: ty = σxy*nx + σyy*ny
    % σyy = λ∂u/∂x + (λ+2μ)∂v/∂y
    % σxy = μ(∂u/∂y + ∂v/∂x)
    
    % 对 u 的贡献（奇数列：1,3,5,...）
    k(2, 1:2:end) = (mu * dphiy * nx + lambda * dphix * ny);
    % 对 v 的贡献（偶数列：2,4,6,...）
    k(2, 2:2:end) = (mu * dphix * nx + (lambda + 2*mu) * dphiy * ny);
end

end
