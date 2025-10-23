function [delta11, delta12, delta21, delta22, xx, yy] = mappings(xcoor, ycoor, xi, eta, type)
%MAPPINGS Isoparametric mapping for several quad variants (serendipity/collapsed)
%
% [delta11, delta12, delta21, delta22, xx, yy] = mappings(xcoor, ycoor, xi, eta, type)
%
% INPUTS
%   xcoor : 1×nen or nen×1   element nodal x-coordinates (nen = #nodes of the chosen type)
%   ycoor : 1×nen or nen×1   element nodal y-coordinates
%   xi    : scalar/array     local ξ  (in [-1,1], unless collapsed mapping specified)
%   eta   : scalar/array     local η  (in [-1,1], unless collapsed mapping specified)
%   type  : integer flag     0  -> Q8 serendipity (8-node)
%                            2  -> collapsed 4-node (corner/edge mapping)
%                            11 -> "Down" collapsed (5-node)
%                            12 -> "Right" collapsed (5-node)
%
% OUTPUTS
%   delta11 = ∂ξ/∂x, delta12 = ∂η/∂x
%   delta21 = ∂ξ/∂y, delta22 = ∂η/∂y
%   xx, yy  = mapped physical coordinates: x(ξ,η), y(ξ,η)
%
% NOTES
% - 本函数返回的是逆雅可比中的条目（∂ξ/∂x 等），便于后续把物理域导数转回参考域/或反之。
% - 为避免在某些退化（collapsed）映射分母为零，按你的原始实现在极限点位做极小量钳制（epss）。
% - 不改变原始形函数与导数公式，只做结构化与注释补强；保持向量化计算。
%
% Author: (W. Huang)
% -------------------------------------------------------------------------

epss = 1.0e-3;  % 极小量，防止退化映射的分母出现 0

% 统一列向量形式，便于矩阵乘
xv = xcoor(:);
yv = ycoor(:);

% ------------------------------
% 依据 type 构造形函数及其一阶导
% ------------------------------
switch type
    case 0
        % ------------------------------
        % Q8 Serendipity (8-node) 元
        % 节点顺序对应你原公式 N1..N8 与 DN* 一致
        % ------------------------------
        N1 = 1/4*(1-xi).*(1-eta).*(-xi-eta-1);
        N2 = 1/4*(1+xi).*(1-eta).*(xi-eta-1);
        N3 = 1/4*(1+xi).*(1+eta).*(xi+eta-1);
        N4 = 1/4*(1-xi).*(1+eta).*(-xi+eta-1);
        N5 = 1/2*(1-xi.^2).*(1-eta);
        N6 = 1/2*(1-eta.^2).*(1+xi);
        N7 = 1/2*(1-xi.^2).*(1+eta);
        N8 = 1/2*(1-eta.^2).*(1-xi);

        DN1X = - (xi/4 - 1/4).*(eta - 1) - ((eta - 1).*(eta + xi + 1))/4;
        DN2X = ((eta - 1).*(eta - xi + 1))/4 - (xi/4 + 1/4).*(eta - 1);
        DN3X = (xi/4 + 1/4).*(eta + 1) + ((eta + 1).*(eta + xi - 1))/4;
        DN4X = (xi/4 - 1/4).*(eta + 1) + ((eta + 1).*(xi - eta + 1))/4;
        DN5X = xi.*(eta - 1);
        DN6X = 1/2 - eta.^2/2;
        DN7X = -xi.*(eta + 1);
        DN8X = eta.^2/2 - 1/2;

        DN1Y = - (xi/4 - 1/4).*(eta - 1) - (xi/4 - 1/4).*(eta + xi + 1);
        DN2Y = (xi/4 + 1/4).*(eta - xi + 1) + (xi/4 + 1/4).*(eta - 1);
        DN3Y = (xi/4 + 1/4).*(eta + 1) + (xi/4 + 1/4).*(eta + xi - 1);
        DN4Y = (xi/4 - 1/4).*(xi - eta + 1) - (xi/4 - 1/4).*(eta + 1);
        DN5Y = xi.^2/2 - 1/2;
        DN6Y = -eta.*(xi + 1);
        DN7Y = 1/2 - xi.^2/2;
        DN8Y = eta.*(xi - 1);

        NN   = [N1(:), N2(:), N3(:), N4(:), N5(:), N6(:), N7(:), N8(:)];
        DNNX = [DN1X(:), DN2X(:), DN3X(:), DN4X(:), DN5X(:), DN6X(:), DN7X(:), DN8X(:)];
        DNNY = [DN1Y(:), DN2Y(:), DN3Y(:), DN4Y(:), DN5Y(:), DN6Y(:), DN7Y(:), DN8Y(:)];

    case 2
        % ------------------------------
        % 右/下方向无限 四边形 (4-node)
        % 注意奇异位置：xi=1 或 eta=-1
        % ------------------------------
        xi2  = xi;  eta2 = eta;
        xi2(xi2==1)     = 1 - epss;
        eta2(eta2==-1)  = -1 + epss;

        N1 = (1-eta2).*(-2*xi2)./((1+eta2).*(1-xi2));
        N2 = (1-eta2).*(1+xi2)./((1+eta2).*(1-xi2));
        N3 = (1+xi2).*(2*eta2)./((1+eta2).*(1-xi2));
        N4 = 4*(-eta2).*xi2./((1+eta2).*(1-xi2));

        DN1X = (2*xi2.*(eta2 - 1))./((eta2 + 1).*(xi2 - 1).^2) - (2*(eta2 - 1))./((eta2 + 1).*(xi2 - 1));
        DN2X = (eta2 - 1)./((eta2 + 1).*(xi2 - 1)) - ((eta2 - 1).*(xi2 + 1))./((eta2 + 1).*(xi2 - 1).^2);
        DN3X = (2*eta2.*(xi2 + 1))./((eta2 + 1).*(xi2 - 1).^2) - (2*eta2)./((eta2 + 1).*(xi2 - 1));
        DN4X = (4*eta2)./((eta2 + 1).*(xi2 - 1)) - (4*eta2.*xi2)./((eta2 + 1).*(xi2 - 1).^2);

        DN1Y = (2*xi2.*(eta2 - 1))./((eta2 + 1).^2.*(xi2 - 1)) - (2*xi2)./((eta2 + 1).*(xi2 - 1));
        DN2Y = (xi2 + 1)./((eta2 + 1).*(xi2 - 1)) - ((eta2 - 1).*(xi2 + 1))./((eta2 + 1).^2.*(xi2 - 1));
        DN3Y = (2*eta2.*(xi2 + 1))./((eta2 + 1).^2.*(xi2 - 1)) - (2*(xi2 + 1))./((eta2 + 1).*(xi2 - 1));
        DN4Y = (4*xi2)./((eta2 + 1).*(xi2 - 1)) - (4*eta2.*xi2)./((eta2 + 1).^2.*(xi2 - 1));

        NN   = [N1(:), N2(:), N3(:), N4(:)];
        DNNX = [DN1X(:), DN2X(:), DN3X(:), DN4X(:)];
        DNNY = [DN1Y(:), DN2Y(:), DN3Y(:), DN4Y(:)];

    case 11
        % ------------------------------
        % "Down" collapsed (5-node)
        % 注意奇异位置：eta = -1
        % ------------------------------
        eta2 = eta;
        eta2(eta2==-1) = -1 + epss;

        N1 = (1-eta2).*(1-xi)./(2*(1+eta2));
        N2 = (1-eta2).*(1+xi)./(2*(1+eta2));
        N3 = (-1+xi+eta2).*(1+xi)./(1+eta2);
        N4 = 2*(1-xi.^2)./(1+eta2);
        N5 = (-1-xi+eta2).*(1-xi)./(1+eta2);

        DN1X = (eta2 - 1)./(2*eta2 + 2);
        DN2X = -(eta2 - 1)./(2*eta2 + 2);
        DN3X = (xi + 1)./(eta2 + 1) + (eta2 + xi - 1)./(eta2 + 1);
        DN5X = (xi - eta2 + 1)./(eta2 + 1) + (xi - 1)./(eta2 + 1);
        DN4X = -(4*xi)./(eta2 + 1);

        DN1Y = (xi - 1)./(2*eta2 + 2) - (2*(eta2 - 1).*(xi - 1))./(2*eta2 + 2).^2;
        DN2Y = (2*(eta2 - 1).*(xi + 1))./(2*eta2 + 2).^2 - (xi + 1)./(2*eta2 + 2);
        DN3Y = (xi + 1)./(eta2 + 1) - ((xi + 1).*(eta2 + xi - 1))./(eta2 + 1).^2;
        DN5Y = - (xi - 1)./(eta2 + 1) - ((xi - 1).*(xi - eta2 + 1))./(eta2 + 1).^2;
        DN4Y = (2*xi.^2 - 2)./(eta2 + 1).^2;

        NN   = [N1(:), N2(:), N3(:), N4(:), N5(:)];
        DNNX = [DN1X(:), DN2X(:), DN3X(:), DN4X(:), DN5X(:)];
        DNNY = [DN1Y(:), DN2Y(:), DN3Y(:), DN4Y(:), DN5Y(:)];

    case 12
        % ------------------------------
        % "Right" collapsed (5-node)
        % 注意奇异位置：xi = 1
        % ------------------------------
        xi2 = xi;
        xi2(xi2==1) = 1 - epss;

        N1 = -(1+xi2+eta).*(1-eta)./(1-xi2);
        N2 = (1+xi2).*(1-eta)./(2*(1-xi2));
        N3 = (1+xi2).*(1+eta)./(2*(1-xi2));
        N4 = (-1+eta-xi2).*(1+eta)./(1-xi2);
        N5 = 2*(1-eta.^2)./(1-xi2);

        DN1X = ((eta - 1).*(eta + xi2 + 1))./(xi2 - 1).^2 - (eta - 1)./(xi2 - 1);
        DN2X = (eta - 1)./(2*xi2 - 2) - (2*(eta - 1).*(xi2 + 1))./(2*xi2 - 2).^2;
        DN3X = (2*(eta + 1).*(xi2 + 1))./(2*xi2 - 2).^2 - (eta + 1)./(2*xi2 - 2);
        DN4X = (eta + 1)./(xi2 - 1) - ((eta + 1).*(xi2 - eta + 1))./(xi2 - 1).^2;
        DN5X = -(2*eta.^2 - 2)./(xi2 - 1).^2;

        DN1Y = - (eta - 1)./(xi2 - 1) - (eta + xi2 + 1)./(xi2 - 1);
        DN2Y = (xi2 + 1)./(2*xi2 - 2);
        DN3Y = -(xi2 + 1)./(2*xi2 - 2);
        DN4Y = (xi2 - eta + 1)./(xi2 - 1) - (eta + 1)./(xi2 - 1);
        DN5Y = (4*eta)./(xi2 - 1);

        NN   = [N1(:), N2(:), N3(:), N4(:), N5(:)];
        DNNX = [DN1X(:), DN2X(:), DN3X(:), DN4X(:), DN5X(:)];
        DNNY = [DN1Y(:), DN2Y(:), DN3Y(:), DN4Y(:), DN5Y(:)];

    otherwise
        error('mappings: unknown type = %d (supported: 0, 2, 11, 12).', type);
end

% --------------- 映射与雅可比 ---------------
% x(ξ,η), y(ξ,η)
xx = NN * xv;
yy = NN * yv;

% 一阶导:
% fxx = ∂x/∂ξ,  fxy = ∂y/∂ξ
% fyx = ∂x/∂η,  fyy = ∂y/∂η
fxx = DNNX * xv;
fxy = DNNX * yv;
fyx = DNNY * xv;
fyy = DNNY * yv;

% 雅可比与其逆
J = fxx .* fyy - fxy .* fyx;

% 告警但不改值，便于定位网格/参数异常
tolJ = 1.0e-12;
if any(abs(J) < tolJ, 'all')
    warning('mappings:JacobianNearlySingular', ...
        'Jacobian magnitude < %.1e detected (%d of %d pts). Check element quality or (xi,eta) range.', ...
         tolJ, nnz(abs(J) < tolJ), numel(J));
end

% 逆雅可比条目
% [∂ξ/∂x  ∂η/∂x;  ∂ξ/∂y  ∂η/∂y] = (1/J) * [ ∂y/∂η  -∂y/∂ξ;  -∂x/∂η  ∂x/∂ξ ]
delta11 =  fyy ./ J;   % ∂ξ/∂x
delta12 = -fxy ./ J;   % ∂η/∂x
delta21 = -fyx ./ J;   % ∂ξ/∂y
delta22 =  fxx ./ J;   % ∂η/∂y
end
