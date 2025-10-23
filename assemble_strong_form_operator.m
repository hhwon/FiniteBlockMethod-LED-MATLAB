function k = assemble_strong_form_operator(phi, dphix, dphiy, dphixx, dphixy, dphiyy, x, ndof, E, nu, plane_strain, axisymmetric, r)
%ASSEMBLE_STRONG_FORM_OPERATOR 装配强形式配点法的刚度矩阵行
%
% k = assemble_strong_form_operator(phi, dphix, dphiy, dphixx, dphixy, dphiyy, x, ndof, E, nu, plane_strain, axisymmetric, r)
%
% INPUTS
%   phi    : 1×n_nodes 行向量，配点处的基函数值
%   dphix  : 1×n_nodes 行向量，配点处 ∂φ/∂x (轴对称时为 ∂φ/∂r)
%   dphiy  : 1×n_nodes 行向量，配点处 ∂φ/∂y (轴对称时为 ∂φ/∂z)
%   dphixx : 1×n_nodes 行向量，配点处 ∂²φ/∂x² (轴对称时为 ∂²φ/∂r²)
%   dphixy : 1×n_nodes 行向量，配点处 ∂²φ/∂x∂y (轴对称时为 ∂²φ/∂r∂z)
%   dphiyy : 1×n_nodes 行向量，配点处 ∂²φ/∂y² (轴对称时为 ∂²φ/∂z²)
%   x      : 标量，配点的 x 坐标（可用于体力项，此处未用）
%   ndof   : 标量，节点总数（注意：全局自由度数 = 2*ndof）
%   E      : 标量，杨氏模量
%   nu     : 标量，泊松比
%   plane_strain : 逻辑值，是否为平面应变（默认 true）
%   axisymmetric : 逻辑值，是否为轴对称问题（默认 false）
%   r      : 标量，配点的径向坐标（仅轴对称时使用，必须 > 0）
%
% OUTPUTS
%   k : 2×(2*ndof) 矩阵，对应该配点的两个方程（x 和 y 方向，或 r 和 z 方向）
%       列按交错顺序排列：[u1, v1, u2, v2, ..., un, vn]
%
% DESCRIPTION
%   本函数实现弹性力学的强形式 Navier 方程：
%   
%   平面应变/应力：
%     (λ+2μ)∂²u/∂x² + μ∂²u/∂y² + (λ+μ)∂²v/∂x∂y = 0
%     (λ+2μ)∂²v/∂y² + μ∂²v/∂x² + (λ+μ)∂²u/∂x∂y = 0
%
%   轴对称（r-z 坐标，u=径向位移，w=轴向位移）：
%     (λ+2μ)∂²u/∂r² + μ∂²u/∂z² + (λ+μ)∂²w/∂r∂z + (λ+2μ)/r·∂u/∂r - (λ+2μ)/r²·u = 0
%     (λ+2μ)∂²w/∂z² + μ∂²w/∂r² + (λ+μ)∂²u/∂r∂z + μ/r·∂w/∂r = 0
%
%   注意：当 r→0 时，使用 L'Hospital 法则：φ/r → dφ/dr (即 dphix)
%
%   其中：λ = Eν/((1+ν)(1-2ν)), μ = E/(2(1+ν)) (平面应变)
%
% See also: assemble_traction_bc_operator
%
% Author: (W. Huang)
% -------------------------------------------------------------------------

if nargin < 11
    plane_strain = true;
end
if nargin < 12
    axisymmetric = false;
end
if nargin < 13
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

% 注意：ndof 是节点数，但全局自由度是 2*ndof（交错排列 u1,v1,u2,v2,...）
n_nodes = ndof;
n_dof_global = 2 * n_nodes;

% 预分配刚度矩阵行 (2行 × 全局自由度数列)
k = zeros(2, n_dof_global);

% 强形式 Navier 方程系数
c1 = lambda + 2 * mu;  % (λ+2μ)
c2 = mu;                % μ
c3 = lambda + mu;       % (λ+μ)

if axisymmetric
    % ========================================================================
    % 轴对称 Navier 方程（r-z 坐标系，u=径向位移，w=轴向位移）
    % ========================================================================
    % r 方向: c1*∂²u/∂r² + c2*∂²u/∂z² + c3*∂²w/∂r∂z + c1/r*∂u/∂r - c1/r²*u = 0
    % z 方向: c1*∂²w/∂z² + c2*∂²w/∂r² + c3*∂²u/∂r∂z + c2/r*∂w/∂r = 0
    %
    % 特殊处理：当 r→0 时，使用 L'Hospital 法则
    %   lim(r→0) φ/r = ∂φ/∂r (即 dphix)
    %   lim(r→0) ∂u/∂r/r = ∂²u/∂r² (即 dphixx)
    % ========================================================================
    
    % 数值稳定性阈值
    r_tol = 1e-5;
    
    if abs(r) < r_tol
        % 当 r ≈ 0 时的极限情况（对称轴上）
        % r 方向: c1*∂²u/∂r² + c2*∂²u/∂z² + c3*∂²w/∂r∂z + c1*∂²u/∂r² - c1*∂u/∂r = 0
        %       = 2*c1*∂²u/∂r² + c2*∂²u/∂z² + c3*∂²w/∂r∂z - c1*∂u/∂r = 0
        % 对 u 的贡献（奇数列）
        k(1, 1:2:end) = 2*c1*dphixx + c2*dphiyy - c1*dphix;
        % 对 w 的贡献（偶数列）
        k(1, 2:2:end) = c3*dphixy;
        
        % z 方向: c1*∂²w/∂z² + c2*∂²w/∂r² + c3*∂²u/∂r∂z + c2*∂²w/∂r² = 0
        %       = c1*∂²w/∂z² + 2*c2*∂²w/∂r² + c3*∂²u/∂r∂z = 0
        % 对 u 的贡献（奇数列）
        k(2, 1:2:end) = c3*dphixy;
        % 对 w 的贡献（偶数列）
        k(2, 2:2:end) = 2*c2*dphixx + c1*dphiyy;
    else
        % 一般情况（r > 0）
        % r 方向方程
        % 对 u 的贡献（奇数列）
        k(1, 1:2:end) = c1*dphixx + c2*dphiyy + c1/r*dphix - c1/(r^2)*phi;
        % 对 w 的贡献（偶数列）
        k(1, 2:2:end) = c3*dphixy;
        
        % z 方向方程
        % 对 u 的贡献（奇数列）
        k(2, 1:2:end) = c3*dphixy;
        % 对 w 的贡献（偶数列）
        k(2, 2:2:end) = c2*dphixx + c1*dphiyy + c2/r*dphix;
    end
else
    % ========================================================================
    % 平面应变/平面应力（x-y 坐标系，u=x方向位移，v=y方向位移）
    % ========================================================================
    % x 方向方程: c1*∂²u/∂x² + c2*∂²u/∂y² + c3*∂²v/∂x∂y = 0
    % 对 u 的贡献（奇数列：1,3,5,...）
    k(1, 1:2:end) = c1 * dphixx + c2 * dphiyy;
    % 对 v 的贡献（偶数列：2,4,6,...）
    k(1, 2:2:end) = c3 * dphixy;
    
    % y 方向方程: c1*∂²v/∂y² + c2*∂²v/∂x² + c3*∂²u/∂x∂y = 0
    % 对 u 的贡献（奇数列：1,3,5,...）
    k(2, 1:2:end) = c3 * dphixy;
    % 对 v 的贡献（偶数列：2,4,6,...）
    k(2, 2:2:end) = c2 * dphixx + c1 * dphiyy;
end

end
