function print_node_summary(xcoor, ycoor, type)
%PRINT_NODE_SUMMARY 打印节点扩展结果的摘要信息
%
% print_node_summary(xcoor, ycoor, type)
%
% INPUTS
%   xcoor : cell 数组，每个单元的 x 坐标
%   ycoor : cell 数组，每个单元的 y 坐标
%   type  : cell 数组，单元类型标识
%
% EXAMPLE
%   [xcoor, ycoor, type] = expand_nodes(CASE.blocks);
%   print_node_summary(xcoor, ycoor, type);
%
% See also: expand_nodes, visualize_nodes
%
% Author: (W. Huang)
% -------------------------------------------------------------------------

nblock = numel(xcoor);

fprintf('\n');
fprintf('========================================\n');
fprintf('       节点生成摘要\n');
fprintf('========================================\n');
fprintf('共处理块数: %d\n', nblock);
fprintf('========================================\n\n');

total_nodes = 0;
for i = 1:nblock
    nnodes = numel(xcoor{i});
    total_nodes = total_nodes + nnodes;
    
    fprintf('Block %d (type=%d): %d 节点\n', i, type{i}, nnodes);
    fprintf('  X 坐标: [%s]\n', num2str(xcoor{i}, '%.4f '));
    fprintf('  Y 坐标: [%s]\n', num2str(ycoor{i}, '%.4f '));
    
    % 计算边界框
    xmin = min(xcoor{i}); xmax = max(xcoor{i});
    ymin = min(ycoor{i}); ymax = max(ycoor{i});
    fprintf('  边界框: X=[%.4f, %.4f], Y=[%.4f, %.4f]\n', xmin, xmax, ymin, ymax);
    fprintf('\n');
end

fprintf('========================================\n');
fprintf('总节点数: %d\n', total_nodes);
fprintf('========================================\n\n');

end
