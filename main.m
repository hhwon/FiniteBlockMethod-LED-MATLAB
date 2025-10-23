%% ==========================================================================
%% 主程序：使用有限/无限块体法求解弹性力学平面问题
%% ==========================================================================
clear; clc; close all;

%% 1. 加载输入文件
CASE = load_json_case('case4_axisym_cylinder.json');
fprintf('加载输入文件: case4_axisym_cylinder.json\n');
fprintf('离散化参数 nt = %d\n\n', CASE.nt);

%% 2. 单元节点扩展：将4节点扩展为 mappings 所需的节点数
fprintf('步骤 1: 单元节点扩展...\n');
[xcoor, ycoor, type] = expand_nodes(CASE.blocks);
nblock = numel(xcoor);
fprintf('  完成！共 %d 个单元\n\n', nblock);

%% 3. 可视化单元节点
fprintf('步骤 2: 可视化单元节点...\n');
options = struct();
options.showNodeLabels  = true;   % 显示节点编号
options.showBlockLabels = true;   % 显示块标签
options.figName         = '单元节点可视化';
visualize_nodes(xcoor, ycoor, type, options);
fprintf('  完成！\n\n');

%% 4. 打印单元节点摘要
print_node_summary(xcoor, ycoor, type);

%% 5. 生成配点（Collocation Points）
fprintf('步骤 3: 生成配点...\n');
[xx, yy, xt, yt, span] = generate_collocation_points(xcoor, ycoor, type, CASE.nt);
fprintf('  完成！每个单元 %d×%d = %d 个配点\n\n', CASE.nt, CASE.nt, CASE.nt^2);

%% 6. 可视化配点分布
fprintf('步骤 4: 可视化配点分布...\n');
opts = struct();
opts.markerSize = 8;
opts.showGrid   = true;   % 显示网格线
opts.figName    = '配点分布';
visualize_collocation_points(xx, yy, opts);
fprintf('  完成！\n\n');

%% 7. 打印配点统计信息
fprintf('========================================\n');
fprintf('       配点生成统计\n');
fprintf('========================================\n');
fprintf('单元数量: %d\n', nblock);
fprintf('每单元配点: %d×%d = %d\n', CASE.nt, CASE.nt, CASE.nt^2);
fprintf('总配点数: %d\n', nblock * CASE.nt^2);
fprintf('参考空间范围: ξ,η ∈ [%.4f, %.4f]\n', min(span), max(span));
fprintf('========================================\n\n');

%% 8. 计算形函数及其导数
fprintf('步骤 5: 计算形函数及其导数...\n');
[Phi, Dphix, Dphiy, Dphixx, Dphixy, Dphiyy, delta11, delta12, delta21, delta22] = ...
    compute_shape_functions(xcoor, ycoor, type, xt, yt, span);
fprintf('  完成！\n\n');

%% 9. 装配全局算子矩阵
fprintf('步骤 6: 装配全局算子矩阵...\n');
[Aphi, ADphix, ADphiy, ADphixx, ADphixy, ADphiyy] = ...
    assemble_global_operators(Phi, Dphix, Dphiy, Dphixx, Dphixy, Dphiyy);
nnodes = CASE.nt^2;  % 每个单元的节点数
fprintf('  完成！全局矩阵尺寸: %d × %d\n', size(Aphi, 1), size(Aphi, 2));
fprintf('  全局自由度数: %d (单元:%d × 节点/单元:%d)\n\n', nblock*nnodes, nblock, nnodes);

%% 10. 打印形函数计算摘要
fprintf('========================================\n');
fprintf('     形函数计算摘要\n');
fprintf('========================================\n');
fprintf('局部算子 (每个单元):\n');
fprintf('  Phi, Dphix, Dphiy: %d × %d\n', nnodes, nnodes);
fprintf('  Dphixx, Dphixy, Dphiyy: %d × %d\n', nnodes, nnodes);
fprintf('\n');
fprintf('全局算子 (所有单元):\n');
fprintf('  Aphi, ADphix, ADphiy: %d × %d\n', size(Aphi, 1), size(Aphi, 2));
fprintf('  ADphixx, ADphixy, ADphiyy: %d × %d\n', size(ADphixx, 1), size(ADphixx, 2));
fprintf('\n');
fprintf('稀疏度统计:\n');
fprintf('  非零元素: %d / %d (%.2f%%)\n', ...
    nnz(Aphi), numel(Aphi), 100*nnz(Aphi)/numel(Aphi));
fprintf('========================================\n\n');

%% 11. 可视化算子稀疏结构
fprintf('步骤 7: 可视化算子稀疏结构...\n');
visualize_operator_sparsity(Aphi, ADphix, ADphiy);
fprintf('  完成！\n\n');

%% 12. 装配全局刚度矩阵和载荷向量
fprintf('步骤 8: 装配全局刚度矩阵...\n');
[K, f] = assemble_global_stiffness(CASE, Aphi, ADphix, ADphiy, ADphixx, ADphixy, ADphiyy, ...
    xx, yy, xt, yt, nblock, nnodes, delta11, delta12, delta21, delta22);
fprintf('  完成！刚度矩阵尺寸: %d × %d\n', size(K, 1), size(K, 2));
fprintf('  载荷向量长度: %d\n\n', length(f));

%% 13. 打印刚度矩阵统计信息
fprintf('========================================\n');
fprintf('     刚度矩阵统计\n');
fprintf('========================================\n');
fprintf('矩阵尺寸: %d × %d\n', size(K, 1), size(K, 2));
fprintf('总自由度数: %d (每节点2个自由度: u, v)\n', size(K, 1));
fprintf('非零元素: %d / %d (%.2f%%)\n', nnz(K), numel(K), 100*nnz(K)/numel(K));
fprintf('载荷向量非零元素: %d\n', nnz(f));
fprintf('条件数估计: %.2e (可能需要预条件)\n', condest(K));
fprintf('========================================\n\n');

%% 14. 求解线性方程组 K*u = f
fprintf('步骤 9: 求解线性方程组...\n');
fprintf('  使用直接求解器...\n');
tic;
u = K \ f;
solve_time = toc;
fprintf('  求解完成！耗时: %.4f 秒\n\n', solve_time);

%% 15. 打印求解结果摘要
fprintf('========================================\n');
fprintf('     求解结果摘要\n');
fprintf('========================================\n');
fprintf('位移向量长度: %d\n', length(u));
fprintf('最大位移: %.6e\n', max(abs(u)));
fprintf('u 方向最大位移: %.6e\n', max(abs(u(1:2:end))));
fprintf('v 方向最大位移: %.6e\n', max(abs(u(2:2:end))));
fprintf('求解时间: %.4f 秒\n', solve_time);
fprintf('========================================\n\n');

%% 16. 可视化求解结果
fprintf('步骤 10: 可视化求解结果...\n');
vis_options = struct();
vis_options.showUndeformed = true;
vis_options.nlevels = 20;  % 等高线层数
fig = visualize_solution(u, xx, yy, nblock, nnodes, vis_options);
fprintf('  完成！\n\n');

fprintf('✅ 所有计算步骤完成！\n\n');


