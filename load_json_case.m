function CASE = load_json_case(json_file)
% 将 schema 的 JSON 转为模块化主程序使用的 CASE 结构
% 用法：CASE = load_json_case('case_linear.json');

S = jsondecode(fileread(json_file));

%% ---------- 基础参数 ----------
CASE.nt  = S.discretization.node;
CASE.tol = 1.0e-8;

%% ---------- 构建几何块（四边形，双线性映射 type=2） ----------
blk_list = S.geometry.blocks;
nb = numel(blk_list);

% 以 id 从 0 起计；MATLAB 内部改为 1 起计的顺序数组
xcoor = cell(1, nb);
ycoor = cell(1, nb);
type  = cell(1, nb);
bbox  = zeros(nb,4); % [xmin xmax ymin ymax]

for k = 1:nb
    b = blk_list(k);
    V = b.vertices;  % 4×2，按 [BL, BR, TR, TL] 或任意顺时针/逆时针给出
    assert((size(V,1)==4 && size(V,2)==2) || ...
           (size(V,1)==8 && size(V,2)==2) || ...
           (size(V,1)==5 && size(V,2)==2), '每个 block 必须给 4/5/8 个顶点 [x,y].');

    % 如果需要，可在这里按多边形方向做一致化；这里默认你提供的是连贯四边形
    xcoor{k} = V(:,1).';  % 行向量 [x1 x2 x3 x4]
    ycoor{k} = V(:,2).';  % 行向量 [y1 y2 y3 y4]
    type{k}  = b.type;    % 块类型
    bbox(k,:) = [min(xcoor{k}), max(xcoor{k}), min(ycoor{k}), max(ycoor{k})];
end

CASE.blocks = struct('xcoor', xcoor, 'ycoor', ycoor, 'type', type);

% 推断整体边块集合（仅在你想用 'LEFT'/'RIGHT' 这类标签时有用）
xmin = min(bbox(:,1)); xmax = max(bbox(:,2));
ymin = min(bbox(:,3)); ymax = max(bbox(:,4));
tolx = 1e-12;
CASE.sets.LEFT   = find(abs(bbox(:,1) - xmin) <= tolx);
CASE.sets.RIGHT  = find(abs(bbox(:,2) - xmax) <= tolx);
CASE.sets.BOTTOM = find(abs(bbox(:,3) - ymin) <= tolx);
CASE.sets.TOP    = find(abs(bbox(:,4) - ymax) <= tolx);

%% ---------- 线弹性材料（支持默认+分块覆盖） ----------
Ecell = cell(1, nb);
vcell = cell(1, nb);

if isfield(S, 'material') && isfield(S.material, 'model') ...
        && strcmpi(S.material.model, 'linear_elastic')
    % 默认
    if isfield(S.material, 'defaults')
        E_def  = S.material.defaults.E;
        nu_def = S.material.defaults.nu;
    else
        E_def  = 1.0;   % 兜底
        nu_def = 0.30;
    end
    for k = 1:nb
        Ecell{k} = E_def;
        vcell{k} = nu_def;
    end
    % 分块覆盖
    if isfield(S.material, 'per_block')
        pb = S.material.per_block;
        ks = fieldnames(pb);
        for i = 1:numel(ks)
            Ecell{i} = pb.(ks{i}).E;
            vcell{i} = pb.(ks{i}).nu;
        end
    end
else
    % 若误传了各向异性刚度（c11/c22...），这里不做解析，直接回退默认各向同性
    warning('material.model 非 linear_elastic，已回退到默认 E,nu=1,0.3');
    for k = 1:nb, Ecell{k}=1.0; vcell{k}=0.30; end
end

CASE.E = Ecell;
CASE.v = vcell;
CASE.plane_strain = true;
if isfield(S.material, 'plane_strain'), CASE.plane_strain = logical(S.material.plane_strain); end

% 轴对称标志
CASE.axisymmetric = false;
if isfield(S.material, 'axisymmetric')
    CASE.axisymmetric = logical(S.material.axisymmetric);
end

%% ---------- 界面（默认为连续 hybrid：一侧力学平衡 + 另一侧位移相等） ----------
CASE.interfaces = struct([]); ci = 0;
if isfield(S,'interfaces')
    for t = 1:numel(S.interfaces)
        I = S.interfaces(t);
        ci = ci + 1;
        CASE.interfaces(ci).type   = 'continuity';
        CASE.interfaces(ci).pair   = [I.master_blk+1, I.slave_blk+1]; % 1 起计
        CASE.interfaces(ci).sides  = { map_face(I.master_face), map_face(I.slave_face) };
        CASE.interfaces(ci).label  = sprintf('blk %d-%d face %s-%s', ...
            I.master_blk, I.slave_blk, upper(I.master_face), upper(I.slave_face));
    end
end

%% ---------- 边界条件（Dirichlet/Neumann） ----------
bc_tmpl = struct('kind','','blocks',[],'side','','comp',[0 0], ...
                 'value',0,'label','');
CASE.bcs = repmat(bc_tmpl, 0, 1);  % 空但字段固定
cb = 0;

if isfield(S,'boundaryConditions')
    for t = 1:numel(S.boundaryConditions)
        B = S.boundaryConditions(t);
        blk  = B.blk  + 1;     % 1 起计
        side = map_face(B.face);
        bcname = lower(B.bc);
        val    = B.value;

        switch bcname
            case 'dirichlet_ux'
                cb = cb + 1;
                CASE.bcs(cb) = make_dirichlet(blk, side, [1 0], val, 'Dirichlet ux');

            case 'dirichlet_uy'
                cb = cb + 1;
                CASE.bcs(cb) = make_dirichlet(blk, side, [0 1], val, 'Dirichlet uy');

            case 'traction_x'
                cb = cb + 1;
                CASE.bcs(cb) = make_neumann(blk, side, [1 0], val, 'Neumann tx');

            case 'traction_y'
                cb = cb + 1;
                CASE.bcs(cb) = make_neumann(blk, side, [0 1], val, 'Neumann ty');

            otherwise
                error('未知 bc 类型：%s', B.bc);
        end
    end
end

%% ---------- 其他元信息 ----------
CASE.solver = S.solver;
CASE.output = S.output;

% ========== 内部小函数 ==========
    function s = map_face(name)
        switch lower(string(name))
            case 'left',   s = 'L';
            case 'right',  s = 'R';
            case 'bottom', s = 'B';
            case 'top',    s = 'T';
            otherwise, error('未知面名：%s', name);
        end
    end

    function bc = make_dirichlet(blks, side, comp, fh, label)
        bc = struct('kind','Dirichlet','blocks',blks,'side',side, ...
                    'comp',comp,'value',fh,'label',label);
    end

    function bc = make_neumann(blks, side, comp, fh, label)
        bc = struct('kind','Neumann','blocks',blks,'side',side, ...
                    'comp',comp,'value',fh,'label',label);
    end
end
