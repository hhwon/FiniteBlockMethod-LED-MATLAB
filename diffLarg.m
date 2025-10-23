function [Phi,Dphix,Dphiy,Dphixx,Dphixy,Dphiyy] = diffLarg(xt, yt, span, delta11, delta12, delta21, delta22)
%DIFFLARG  基于一维 Lagrange 基函数（x、y 方向张量积）的二维基与导数算子装配
%
% 目的：
%   在给定配点 (xt, yt) 与一维结点 span 的情况下，构造二维张量积基函数
%   φ(x,y)=Lx(x)*Ly(y) 在每个配点的函数值与一阶导数（物理坐标系下），并通过
%   一阶导数算子复合得到二阶导数算子。
%
% 输入参数：
%   xt, yt   : 一维配点（向量）。函数内部默认使用 numel(xt) 作为矩阵阶数，
%              因而通常要求 numel(xt)==numel(yt)==numel(span)。
%   span     : 一维插值结点（向量），用于构造 Lagrange 基（与 larg 的第二参一致）。
%   delta11, delta12, delta21, delta22 :
%              逆雅可比的四个分量，长度应为 numel(xt)，在每个配点 i 处
%              使用：
%                  delta11=∂ξ/∂x, delta12=∂η/∂x,
%                  delta21=∂ξ/∂y, delta22=∂η/∂y。
%
% 输出参数（均为 n×n 方阵，n=numel(xt)）：
%   Phi      : 基函数值矩阵。第 i 行对应第 i 个配点处，对“所有自由度/结点”的基值。
%   Dphix    : 一阶导数算子（物理 x 方向）。左乘自由度向量即得配点处 ∂φ/∂x。
%   Dphiy    : 一阶导数算子（物理 y 方向）。左乘自由度向量即得配点处 ∂φ/∂y。
%   Dphixx   : 二阶导数算子（物理 x 方向），由 Dphix*Dphix 复合近似获得。
%   Dphixy   : 混合二阶导数算子（对称化组合 (Dphix*Dphiy + Dphiy*Dphix)/2）。
%   Dphiyy   : 二阶导数算子（物理 y 方向），由 Dphiy*Dphiy 复合近似获得。
%
% 计算要点（与 larg 的约定一致）：
%   [ux, uxdx] = larg(xt, span)  返回：
%       ux(i,j)   = Lx_i(xt_j)         —— x 方向基函数值
%       uxdx(i,j) = dLx_i/dx |xt_j     —— x 方向一维导数
%   [uy, uydy] = larg(yt, span)  返回：
%       uy(i,j)   = Ly_i(yt_j)         —— y 方向基函数值
%       uydy(i,j) = dLy_i/dy |yt_j     —— y 方向一维导数
%
%   在配点 i 处，二维张量积形式：
%       φ =  Lx .* Ly
%       ∂φ/∂ξ = (dLx/dx_1D) .* Ly
%       ∂φ/∂η =  Lx .* (dLy/dy_1D)
%   再用链式法则从 (ξ,η) 到 (x,y)：
%       ∂φ/∂x = (∂φ/∂ξ)*delta11(i) + (∂φ/∂η)*delta12(i)
%       ∂φ/∂y = (∂φ/∂ξ)*delta21(i) + (∂φ/∂η)*delta22(i)
%
% 说明：
%   1) 这里二阶导数通过一阶算子复合得到（谱/强形式中常用且简洁）。
%   2) 若几何或系数在点间变化剧烈，直接构造二阶算子可能更精确。
%   3) 维度匹配：本函数假设 span 的长度与自由度数一致，且与 xt/yt 保持一致。
%      若在局部出现维度不符错误，请检查 ux/uy 的转置使用是否与期望一致。
% -------------------------------------------------------------------------

Dphix = zeros(numel(xt), numel(xt));   % 预分配：x 方向一阶导算子矩阵（n×n）
Dphiy = zeros(numel(xt), numel(xt));   % 预分配：y 方向一阶导算子矩阵（n×n）
Phi   = zeros(numel(xt), numel(xt));   % 预分配：基函数值矩阵（n×n）

% 基于一维结点 span，分别在 xt、yt 上构造 Lagrange 基与一维导数
[ux, uxdx] = larg(xt, span);           % ux:  Lx(xt)，   uxdx: dLx/dx(xt)
[uy, uydy] = larg(yt, span);           % uy:  Ly(yt)，   uydy: dLy/dy(yt)

% 逐配点装配行算子（第 i 行对应第 i 个配点）
for i = 1:numel(xt)
    % --- 参考系下的导数/基值（张量积） ---
    % udx：∂φ/∂ξ = (dLx/dx_1D).*Ly
    udx = uxdx(:,i)'.*uy(:,i);
    % udy：∂φ/∂η = Lx.*(dLy/dy_1D)
    udy = ux(:,i)'.*uydy(:,i);
    % phi：φ = Lx.*Ly
    phi = ux(:,i)'.*uy(:,i);

    % --- 链式法则：映射到物理坐标系下的一阶导数 ---
    % dphix = ∂φ/∂x
    dphix = delta11(i).*udx(:)+delta12(i)*udy(:);
    % dphiy = ∂φ/∂y
    dphiy = delta21(i).*udx(:)+delta22(i)*udy(:);

    % 将行向量装入对应的第 i 行
    Phi(i,:)=phi(:);
    Dphix(i, :) = dphix;
    Dphiy(i, :) = dphiy;
end

% --- 由一阶算子复合得到二阶算子 ---
Dphixx = Dphix*Dphix;                              % 近似 ∂^2/∂x^2
Dphixy = (Dphix*Dphiy+Dphiy*Dphix)/2;              % 对称化的 ∂^2/(∂x∂y)
Dphiyy = Dphiy*Dphiy;                              % 近似 ∂^2/∂y^2
end
