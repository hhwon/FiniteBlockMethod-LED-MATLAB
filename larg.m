function [u, du] = larg(x, xm)
%LARG  一维 Lagrange 基函数及其一阶导数（重心/Barycentric 公式）
%
% 用途：
%   在给定一维结点 xm 与评估点 x 的情况下，利用重心（barycentric）插值公式
%   构造 Lagrange 基函数值矩阵 u 及其一阶导数矩阵 du。
%   本函数常与 diffLarg 联用以装配 2D 张量积基与导数算子。
%
% 输入参数：
%   x   : 评估点向量（行/列均可），长度为 xn
%   xm  : 插值结点向量（行/列均可），长度为 xmn（也是基函数个数）
%
% 输出参数（均为 xmn×xn 方阵/矩阵）：
%   u   : 基函数值矩阵，u(i,j) = L_i( x(j) )
%   du  : 基函数一阶导矩阵，du(i,j) = dL_i/dx |_{x = x(j)}
%
% 计算要点（重心插值）：
%   - 重心权重：w_i = 1 / Π_{k≠i} (xm_i - xm_k)
%   - 对任意 x：
%       Z_i(x) = w_i / (x - xm_i),   R(x) = Σ_i Z_i(x),   S(x) = Σ_i w_i / (x - xm_i)^2
%       L_i(x) =  Z_i(x) / R(x)
%       dL_i/dx = [ Z_i(x) * S(x) - w_i/(x - xm_i)^2 * R(x) ] / R(x)^2
%   - 命中节点（x ≈ xm_i）时的解析处理（避免除零、保证卡氏性）：
%       L_i(xm_i)=1, L_{k≠i}(xm_i)=0
%       dL_{k≠i}/dx|_{xm_i} = w_k / [ w_i * (xm_i - xm_k) ],
%       dL_i/dx|_{xm_i}     = -Σ_{k≠i} dL_{k≠i}/dx
%
% 说明：
%   1) 形状规范：输出 u、du 的尺寸均为 xmn×xn，以“基函数索引为行、评估点为列”的约定与 diffLarg 对齐。
%   2) 数值健壮性：对 x≈xm_i 采用自适应容差识别“命中”，并用闭式导数替代通式，避免数值不稳。
%   3) 复杂度：权重 w 的预处理 O(xmn^2)（一次性）；求值 O(xmn·xn)。
%   4) 约束：若 xm 含重复结点会导致权重发散，函数将报错。
% -------------------------------------------------------------------------

    % 统一形状：xm 作为列向量、x 作为行向量，便于隐式扩展
    x  = x(:).';   % 1 × xn
    xm = xm(:);    % xmn × 1

    xmn = numel(xm);  % 节点数（基函数个数）
    xn  = numel(x);   % 评估点数

    % --- 基本健壮性检查：重复节点将导致权重发散 ---
    if any(diff(sort(xm)) == 0)
        error('larg:duplicateNodes', 'xm contains duplicate nodes.');
    end

    % --- 计算重心权重 w_i = 1 / Π_{k≠i} (xm_i - xm_k) ---
    Xdiff = xm - xm.';                  % xmn × xmn, 第(i,j)=xm_i - xm_j
    Xdiff(1:xmn+1:end) = 1;             % 对角置 1，避免乘积为 0
    w = 1 ./ prod(Xdiff, 2);            % xmn × 1

    % --- 与评估点的差分矩阵 Δ = x - xm ---
    % 隐式扩展： (xmn×1) 对 (1×xn) 扩展为 xmn×xn
    Delta = x - xm;                     % xmn × xn

    % --- 识别“命中节点”的列（x(j)=xm(i)）并用有限容差处理 ---
    xm_sorted = sort(xm);
    if xmn > 1
        minsep = min(abs(diff(xm_sorted)));
    else
        minsep = 1;
    end
    tol = eps(max(1, minsep)) * 16;     % 可按需调整

    hitMask = abs(Delta) <= tol;        % xmn × xn
    anyHit  = any(hitMask, 1);          % 1 × xn，哪些列命中节点

    % --- 通用列（未命中节点）采用重心公式 ---
    % Z = w./(x - xm),  R = Σ Z,  S = Σ w./(x - xm).^2
    Z = w ./ Delta;                     % xmn × xn
    R = sum(Z, 1);                      % 1 × xn
    S = sum(w ./ (Delta.^2), 1);        % 1 × xn

    % L_i(x) 与导数（先按一般情形计算）
    u  = Z ./ R;                                                % xmn × xn
    du = (Z .* S - (w ./ (Delta.^2)) .* R) ./ (R.^2);           % xmn × xn

    % --- 处理命中列：强制插值“卡氏性”与闭式导数 ---
    if any(anyHit)
        hitCols = find(anyHit);          % 这些列有 x(j)≈某个 xm(i0)
        for c = hitCols
            rows = find(hitMask(:, c));  % 通常只会有一个 i0
            i0 = rows(1);

            % 基函数值：严格的 Kronecker delta
            u(:, c) = 0;
            u(i0, c) = 1;

            % 基函数导数（闭式表达）
            % 对 k≠i0:  du(k,c) = w(k) / ( w(i0) * (xm(i0) - xm(k)) )
            du(:, c) = 0;
            idx = true(xmn,1); idx(i0) = false;
            du(idx, c) = w(idx) ./ ( w(i0) * (xm(i0) - xm(idx)) );

            % 对 i0: du(i0,c) = -sum_{k≠i0} du(k,c)
            du(i0, c) = -sum(du(idx, c));
        end
    end

    % --- 与 diffLarg 约定一致的返回形状 ---
    % u、du 已为 xmn×xn，无需 reshape

end
