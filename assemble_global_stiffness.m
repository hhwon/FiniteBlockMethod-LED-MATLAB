function [K, f] = assemble_global_stiffness(CASE, Aphi, ADphix, ADphiy, ADphixx, ADphixy, ADphiyy, xx, yy, xt, yt, nblock, nnodes, delta11, delta12, delta21, delta22)
%ASSEMBLE_GLOBAL_STIFFNESS 装配全局刚度矩阵和载荷向量
%
% [K, f] = assemble_global_stiffness(CASE, Aphi, ADphix, ADphiy, ADphixx, ADphixy, ADphiyy, xx, yy, xt, yt, nblock, nnodes, delta11, delta12, delta21, delta22)
%
% INPUTS
%   CASE    : 算例结构体，包含材料参数、边界条件等
%   Aphi    : 全局基函数值矩阵
%   ADphix  : 全局 x 方向一阶导数算子
%   ADphiy  : 全局 y 方向一阶导数算子
%   ADphixx : 全局 x 方向二阶导数算子
%   ADphixy : 全局混合二阶导数算子
%   ADphiyy : 全局 y 方向二阶导数算子
%   xx      : cell 数组，每个单元的配点 x 坐标
%   yy      : cell 数组，每个单元的配点 y 坐标
%   xt      : nt×nt 矩阵，参考空间 ξ 坐标
%   yt      : nt×nt 矩阵，参考空间 η 坐标
%   nblock  : 单元数量
%   nnodes  : 每个单元的配点数
%   delta11 : cell 数组，逆雅可比分量 ∂ξ/∂x
%   delta12 : cell 数组，逆雅可比分量 ∂η/∂x
%   delta21 : cell 数组，逆雅可比分量 ∂ξ/∂y
%   delta22 : cell 数组，逆雅可比分量 ∂η/∂y
%
% OUTPUTS
%   K : (2*nblock*nnodes) × (2*nblock*nnodes) 全局刚度矩阵
%   f : (2*nblock*nnodes) × 1 全局载荷向量
%
% DESCRIPTION
%   本函数装配全局刚度矩阵，包括：
%   1. 内部配点：强形式 Navier 方程
%   2. 边界配点：Dirichlet 或 Neumann 边界条件
%   3. 界面配点：连续性条件（位移相等 + 应力平衡）
%
%   自由度排列：交错排列 [u1, v1, u2, v2, ..., un, vn]
%
% See also: assemble_strong_form_operator, assemble_traction_bc_operator, compute_boundary_normal
%
% Author: (W. Huang)
% -------------------------------------------------------------------------

fprintf('  开始装配全局刚度矩阵...\n');

% 初始化全局刚度矩阵和载荷向量
n_nodes_total = nblock * nnodes;   % 总节点数
ndof_total = 2 * n_nodes_total;    % 总自由度数 (u 和 v)
K = zeros(ndof_total, ndof_total);
f = zeros(ndof_total, 1);

%% 遍历所有单元和配点
for iblock = 1:nblock
    fprintf('    处理单元 %d/%d...\n', iblock, nblock);
    
    % 当前单元的材料参数
    E = CASE.E{iblock};
    nu = CASE.v{iblock};
    
    % 当前单元的逆雅可比矩阵（用于计算边界法向量）
    delta11_block = delta11{iblock};
    delta12_block = delta12{iblock};
    delta21_block = delta21{iblock};
    delta22_block = delta22{iblock};
    
    for inode = 1:nnodes
        % 全局配点编号
        num = (iblock - 1) * nnodes + inode;
        
        % 当前配点在参考空间的位置
        xi = xt(inode);
        eta = yt(inode);
        
        % 提取当前配点的算子行
        phi_row    = Aphi(num, :);
        dphix_row  = ADphix(num, :);
        dphiy_row  = ADphiy(num, :);
        dphixx_row = ADphixx(num, :);
        dphixy_row = ADphixy(num, :);
        dphiyy_row = ADphiyy(num, :);
        
        % 当前配点的物理坐标（从 cell 数组中获取当前单元的坐标矩阵）
        xx_block = xx{iblock};  % nt×nt 矩阵
        yy_block = yy{iblock};  % nt×nt 矩阵
        x_pos = xx_block(inode);  % 线性索引
        y_pos = yy_block(inode);  % 线性索引
        
        % 标记是否已处理（边界或界面）
        is_processed = false;
        
        % 边界判断：配点在参考空间边界上
        tol = 1e-5;
        
        % 判断配点是否在某个边界上（这些变量在后续的界面和边界检查中都会用到）
        on_left   = abs(xi - (-1)) < tol;
        on_right  = abs(xi - 1) < tol;
        on_bottom = abs(eta - (-1)) < tol;
        on_top    = abs(eta - 1) < tol;
        on_any_boundary = on_left || on_right || on_bottom || on_top;
        
        %% 1. 检查是否为界面配点
        % ====================================================================
        % 界面连续性条件处理（修正版）
        % ====================================================================
        % 在界面上需要满足：
        %   1. 位移连续：u_master - u_slave = 0, v_master - v_slave = 0
        %   2. 应力平衡：σ_master·n_master + σ_slave·n_slave = 0
        %      （即：t_master + t_slave = 0，牵引力大小相等方向相反）
        %
        % 实现方式（关键修正）：
        %   - 在主侧配点：施加应力平衡条件
        %     方程：σ_master·n_master + σ_slave·n_slave = 0
        %     需要同时包含主侧和从侧的应力算子
        %
        %   - 在从侧配点：施加位移连续条件
        %     方程：u_slave - u_master = 0, v_slave - v_master = 0
        %     需要同时包含从侧和主侧的位移
        %
        % 假设：界面上主从节点一一对应（网格对齐）
        % ====================================================================
        
        % 检查当前配点是否在界面上
        is_on_interface = false;
        interface_master = false;  % 是否为主侧
        interface_info = struct();  % 存储界面信息
        
        if isfield(CASE, 'interfaces') && ~isempty(CASE.interfaces)
            for iface = 1:numel(CASE.interfaces)
                ifc = CASE.interfaces(iface);  % 结构体数组，使用圆括号
                
                % 获取主从单元编号（已经是 1-based）
                master_blk = ifc.pair(1);
                slave_blk = ifc.pair(2);
                master_face = ifc.sides{1};  % 'L', 'R', 'B', 'T'
                slave_face = ifc.sides{2};
                
                % 检查是否为主侧界面配点
                if master_blk == iblock
                    % 检查是否在主侧界面上
                    if (strcmp(master_face, 'L') && on_left) || ...
                       (strcmp(master_face, 'R') && on_right) || ...
                       (strcmp(master_face, 'B') && on_bottom) || ...
                       (strcmp(master_face, 'T') && on_top)
                        is_on_interface = true;
                        interface_master = true;
                        
                        % 保存界面信息
                        interface_info.master_blk = master_blk;
                        interface_info.slave_blk = slave_blk;
                        interface_info.master_face = master_face;
                        interface_info.slave_face = slave_face;
                        interface_info.inode_local = inode;  % 当前节点在本单元的局部编号
                        
                        % 动态计算主侧法向量（使用已有的雅可比矩阵信息）
                        delta11_master = delta11{master_blk};
                        delta12_master = delta12{master_blk};
                        delta21_master = delta21{master_blk};
                        delta22_master = delta22{master_blk};
                        [nx_m, ny_m] = compute_boundary_normal(xi, eta, delta11_master, delta12_master, ...
                                                                   delta21_master, delta22_master, xt, yt, master_face);
                        interface_info.n_master = [nx_m, ny_m];
                        
                        % 从侧法向量（指向主侧，与主侧相反）
                        interface_info.n_slave = -interface_info.n_master;
                        
                        break;
                    end
                end
                
                % 检查是否为从侧界面配点
                if slave_blk == iblock
                    % 检查是否在从侧界面上
                    if (strcmp(slave_face, 'L') && on_left) || ...
                       (strcmp(slave_face, 'R') && on_right) || ...
                       (strcmp(slave_face, 'B') && on_bottom) || ...
                       (strcmp(slave_face, 'T') && on_top)
                        is_on_interface = true;
                        interface_master = false;
                        
                        % 保存界面信息
                        interface_info.master_blk = master_blk;
                        interface_info.slave_blk = slave_blk;
                        interface_info.master_face = master_face;
                        interface_info.slave_face = slave_face;
                        interface_info.inode_local = inode;  % 当前节点在本单元的局部编号
                        
                        % 动态计算主侧法向量（使用已有的雅可比矩阵信息）
                        % 角点情况：需要两个法向量
                        delta11_master = delta11{master_blk};
                        delta12_master = delta12{master_blk};
                        delta21_master = delta21{master_blk};
                        delta22_master = delta22{master_blk};
                        
                        % 计算每个边界面的法向量
                        [nx_m1, ny_m1] = compute_boundary_normal(xi, eta, delta11_master, delta12_master, ...
                                                                     delta21_master, delta22_master, xt, yt, master_face);
                        [nx_m2, ny_m2] = compute_boundary_normal(xi, eta, delta11_master, delta12_master, ...
                                                                     delta21_master, delta22_master, xt, yt, slave_face);
                        
                        % 存储两个法向量（每行一个）
                        interface_info.n_master = [nx_m1, ny_m1; nx_m2, ny_m2];
                        
                        % 从侧法向量（与主侧相反）
                        interface_info.n_slave = -interface_info.n_master;
                        
                        break;
                    end
                end
            end
        end
        
        % 如果在界面上，处理界面条件
        if is_on_interface
            
            
            if interface_master
                % 计算主从节点的全局编号（假设界面对齐，节点一一对应）
                slave_blk_inode_local = find(abs(xx{interface_info.slave_blk}-x_pos)<=tol & abs(yy{interface_info.slave_blk}-y_pos)<=tol);
                num_slave = (interface_info.slave_blk - 1) * nnodes + slave_blk_inode_local;
                num_master = (interface_info.master_blk - 1) * nnodes + interface_info.inode_local;
                % ============================================================
                % 主侧配点：施加应力平衡条件
                % ============================================================
                % 方程：t_master + t_slave = 0
                %       σ_master·n_master + σ_slave·n_slave = 0
                %
                % 需要同时包含主侧和从侧的应力算子
                % ============================================================
                
                % 主侧材料参数和算子
                E_master = CASE.E{interface_info.master_blk};
                nu_master = CASE.v{interface_info.master_blk};
                
                % 提取主侧配点（当前配点）的形函数
                phi_master = Aphi(num_master, :);
                dphix_master = ADphix(num_master, :);
                dphiy_master = ADphiy(num_master, :);
                
                % 主侧物理坐标
                xx_master = xx{interface_info.master_blk};
                yy_master = yy{interface_info.master_blk};
                x_master = xx_master(interface_info.inode_local);
                
                % 主侧应力算子：σ_master·n_master
                k_master = assemble_traction_bc_operator(phi_master, dphix_master, dphiy_master, ...
                    x_master, n_nodes_total, interface_info.n_master, E_master, nu_master, ...
                    CASE.plane_strain, CASE.axisymmetric, x_master);
                
                % 从侧材料参数和算子
                E_slave = CASE.E{interface_info.slave_blk};
                nu_slave = CASE.v{interface_info.slave_blk};
                
                % 提取从侧对应配点的形函数
                phi_slave = Aphi(num_slave, :);
                dphix_slave = ADphix(num_slave, :);
                dphiy_slave = ADphiy(num_slave, :);
                
                % 从侧物理坐标（应该与主侧相同）
                xx_slave = xx{interface_info.slave_blk};
                yy_slave = yy{interface_info.slave_blk};
                x_slave = xx_slave(slave_blk_inode_local);
                
                % 从侧应力算子：σ_slave·n_slave
                k_slave = assemble_traction_bc_operator(phi_slave, dphix_slave, dphiy_slave, ...
                    x_slave, n_nodes_total, interface_info.n_slave, E_slave, nu_slave, ...
                    CASE.plane_strain, CASE.axisymmetric, x_slave);
                
                % 应力平衡方程：t_master + t_slave = 0
                % 即：σ_master·n_master + σ_slave·n_slave = 0
                K(2*num-1:2*num, :) = k_master + k_slave;
                f(2*num-1:2*num) = [0; 0];  % 右端项为 0（应力平衡）
                
            else
                % 计算主从节点的全局编号（假设界面对齐，节点一一对应）
                master_blk_inode_local = find(abs(xx{interface_info.master_blk}-x_pos)<=tol & abs(yy{interface_info.master_blk}-y_pos)<=tol);
                num_master = (interface_info.master_blk - 1) * nnodes + master_blk_inode_local;
                num_slave = (interface_info.slave_blk - 1) * nnodes + interface_info.inode_local;
                % ============================================================
                % 从侧配点：施加位移连续条件
                % ============================================================
                % 方程：u_slave - u_master = 0
                %       v_slave - v_master = 0
                %
                % 这需要同时包含从侧和主侧的位移
                % ============================================================
                
                % u 方向：u_slave - u_master = 0
                K(2*num-1, :) = 0;  % 清空该行
                K(2*num-1, 2*num_slave-1) = 1.0;   % u_slave 的系数
                K(2*num-1, 2*num_master-1) = -1.0; % -u_master 的系数
                
                % v 方向：v_slave - v_master = 0
                K(2*num, :) = 0;  % 清空该行
                K(2*num, 2*num_slave) = 1.0;       % v_slave 的系数
                K(2*num, 2*num_master) = -1.0;     % -v_master 的系数
                
                f(2*num-1:2*num) = [0; 0];  % 右端项为 0（位移连续）
            end
            is_processed = true;
        end
        
        if is_processed
            continue;
        end
        
        %% 2. 检查是否为边界配点
        % ====================================================================
        % 关键原理：强形式配点法的边界条件处理
        % ====================================================================
        % 强形式配点法要求：
        %   - 内部配点：满足强形式控制方程（Navier方程）
        %   - 边界配点：满足边界条件（位移或应力）
        %
        % 对于边界配点，有三种情况：
        %   1. 显式 Dirichlet 边界条件（指定位移）
        %   2. 显式 Neumann 边界条件（指定应力/牵引力）
        %   3. 自由边界（未显式指定 = 隐式零应力边界）
        %
        % 重要：所有边界配点都必须施加边界条件，包括自由边界！
        %       自由边界 = 应力为零的Neumann边界条件
        % ====================================================================
        
        % 使用前面定义的 on_left, on_right, on_bottom, on_top, on_any_boundary
        
        if ~on_any_boundary
            % 不在边界上，继续处理内部点或界面点
        else
            % 在边界上，需要应用边界条件（显式或默认）
            
            % 初始化边界条件
            k_bc = zeros(2, ndof_total);
            f_bc = zeros(2, 1);
            k_dirichlet = zeros(2, ndof_total);
            has_prescribed_bc = false;  % 是否有显式指定的边界条件
            
            % 记录当前配点所在的所有边界
            boundaries = {};
            if on_left,   boundaries{end+1} = 'L'; end
            if on_right,  boundaries{end+1} = 'R'; end
            if on_bottom, boundaries{end+1} = 'B'; end
            if on_top,    boundaries{end+1} = 'T'; end
            
            % 获取当前单元的物理坐标网格
            xx_block = xx{iblock};
            yy_block = yy{iblock};

            % 角点：计算所有边界的法向量并取平均
            n_avg = zeros(1, 2);
            for ib = 1:numel(boundaries)
                face = boundaries{ib};
                [nx, ny] = compute_boundary_normal(xi, eta, delta11_block, delta12_block, ...
                                                       delta21_block, delta22_block, xt, yt, face);
                n_avg = n_avg + [nx, ny];
            end
            % 归一化
            n_norm = norm(n_avg);
            if n_norm > eps
                normal = n_avg / n_norm;
            else
                % 如果法向量和为零（不应该发生），使用第一个边界
                face = boundaries{1};
                [nx, ny] = compute_boundary_normal(xi, eta, delta11_block, delta12_block, ...
                                                       delta21_block, delta22_block, xt, yt, face);
                normal = [nx, ny];
                warning('角点法向量和为零，使用第一个边界法向量');
            end
            k_traction = assemble_traction_bc_operator(phi_row, dphix_row, dphiy_row, ...
                    x_pos, n_nodes_total, normal, E, nu, CASE.plane_strain, ...
                    CASE.axisymmetric, x_pos);
            k_bc = k_traction;
            
            % 遍历所有边界条件，查找是否有显式指定的边界条件
            for ibc = 1:numel(CASE.bcs)

                bc = CASE.bcs(ibc);
                
                % 检查是否应用到当前单元
                if ~ismember(iblock, bc.blocks)
                    continue;
                end
                
                % 检查是否在指定边界上
                if ~any(strcmp(bc.side, boundaries))
                    continue;
                end
                
                % 动态计算边界法向量（使用已有的雅可比矩阵信息）
                idx = find(strcmp(bc.side, boundaries), 1);
                face = boundaries{idx};
                [nx, ny] = compute_boundary_normal(xi, eta, delta11_block, delta12_block, ...
                                                       delta21_block, delta22_block, xt, yt, face);
                normal = [nx, ny];

                if strcmp(bc.kind, 'Dirichlet')
                    
                    
                    % u 方向
                    if bc.comp(1) == 1
                        k_bc(1, :) = 0;
                        k_bc(1, 2*num-1:2*num) = [1, 0];
                        f_bc(1) = bc.value;
                    end
                    
                    % v 方向
                    if bc.comp(2) == 1
                        k_bc(2, :) = 0;
                        k_bc(2, 2*num-1:2*num) = [0, 1];
                        f_bc(2) = bc.value;
                    end
                    
                elseif strcmp(bc.kind, 'Neumann')
                    % Neumann 边界条件：指定牵引力
                    k_traction = assemble_traction_bc_operator(phi_row, dphix_row, dphiy_row, ...
                        x_pos, n_nodes_total, normal, E, nu, CASE.plane_strain, ...
                        CASE.axisymmetric, x_pos);
                    
                    % 载荷向量
                    if bc.comp(1) == 1
                        k_bc(1, :) = k_traction(1, :);
                        f_bc(1) = bc.value;
                    end
                    if bc.comp(2) == 1
                        k_bc(2, :) = k_traction(2, :);
                        f_bc(2) = bc.value;
                    end
                end
            end

            % 应用边界条件
            K(2*num-1:2*num, :) = k_bc;
            f(2*num-1:2*num) = f_bc;
            is_processed = true;
        end
        
        if is_processed
            continue;
        end
        
        %% 3. 内部配点：应用强形式 Navier 方程
        k_strong = assemble_strong_form_operator(phi_row, dphix_row, dphiy_row, ...
            dphixx_row, dphixy_row, dphiyy_row, x_pos, n_nodes_total, E, nu, ...
            CASE.plane_strain, CASE.axisymmetric, x_pos);
        
        K(2*num-1:2*num, :) = k_strong;
        f(2*num-1:2*num) = [0; 0];  % 无体力
    end
end

fprintf('  刚度矩阵装配完成！\n');

end
