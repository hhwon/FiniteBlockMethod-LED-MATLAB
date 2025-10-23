function [xcoor_expanded, ycoor_expanded, type_expanded] = expand_nodes(blocks)
%EXPAND_NODES 将输入的4节点单元扩展为 mappings 所需的节点数
%
% [xcoor_expanded, ycoor_expanded, type_expanded] = expand_nodes(blocks)
%
% INPUTS
%   blocks : 结构体，包含字段 xcoor, ycoor, type (均为 cell 数组)
%            xcoor{i}, ycoor{i} : 输入的4个角点坐标 (行向量)
%            type{i}            : 单元类型标识
%
% OUTPUTS
%   xcoor_expanded : cell 数组，扩展后的 x 坐标
%   ycoor_expanded : cell 数组，扩展后的 y 坐标
%   type_expanded  : cell 数组，单元类型（保持不变）
%
% ELEMENT TYPES
%   type=0  : Q8 Serendipity (8-node) - 4角点 -> 4角点 + 4边中点
%   type=2  : 4-node 无限单元 - 保持4节点
%   type=11 : "Down" collapsed (5-node) - 4角点 -> 5节点
%   type=12 : "Right" collapsed (5-node) - 4角点 -> 5节点
%
% EXAMPLE
%   CASE = load_json_case('case_linear.json');
%   [xcoor, ycoor, type] = expand_nodes(CASE.blocks);
%
% See also: mappings, load_json_case
%
% Author: (W. Huang)
% -------------------------------------------------------------------------

nblock = numel(blocks);

xcoor_expanded = cell(1, nblock);
ycoor_expanded = cell(1, nblock);
type_expanded  = cell(1, nblock);

for i = 1:nblock
    % 从输入获取节点坐标
    xcoor_input = blocks(i).xcoor;  % 行向量
    ycoor_input = blocks(i).ycoor;  % 行向量
    type_input  = blocks(i).type;
    
    % 获取输入节点数
    n_input = length(xcoor_input);
    
    % 确定该单元类型所需的节点数
    switch type_input
        case 0
            n_required = 8;  % Q8 Serendipity
        case 2
            n_required = 4;  % 4-node 无限单元
        case 11
            n_required = 5;  % "Down" collapsed
        case 12
            n_required = 5;  % "Right" collapsed
        otherwise
            error('expand_nodes:UnknownType', '未知的单元类型: type=%d', type_input);
    end
    
    % 判断是否需要扩展
    if n_input == n_required
        % 输入节点数已满足要求，直接使用（支持曲边单元）
        xcoor_expanded{i} = xcoor_input;
        ycoor_expanded{i} = ycoor_input;
        fprintf('  Block %d (type=%d): 输入节点数=%d，已满足要求，无需扩展\n', ...
            i, type_input, n_input);
    elseif n_input == 4
        % 输入4个角点，需要扩展
        switch type_input
            case 0
                % Q8 Serendipity (8-node)
                [xcoor_expanded{i}, ycoor_expanded{i}] = expand_to_8nodes(xcoor_input, ycoor_input);
                
            case 2
                % 4-node 无限单元，保持4节点
                xcoor_expanded{i} = xcoor_input;
                ycoor_expanded{i} = ycoor_input;
                
            case 11
                % "Down" collapsed (5-node)
                [xcoor_expanded{i}, ycoor_expanded{i}] = expand_to_5nodes_down(xcoor_input, ycoor_input);
                
            case 12
                % "Right" collapsed (5-node)
                [xcoor_expanded{i}, ycoor_expanded{i}] = expand_to_5nodes_right(xcoor_input, ycoor_input);
        end
        fprintf('  Block %d (type=%d): 从%d个节点扩展到%d个节点\n', ...
            i, type_input, n_input, length(xcoor_expanded{i}));
    else
        % 节点数不匹配
        error('expand_nodes:NodeCountMismatch', ...
            'Block %d (type=%d): 输入节点数=%d，应为4（需扩展）或%d（已完整）', ...
            i, type_input, n_input, n_required);
    end
    
    type_expanded{i} = type_input;
end

end

%% ========================================================================
%% 子函数：具体的节点扩展逻辑
%% ========================================================================

function [xout, yout] = expand_to_8nodes(xin, yin)
% 将4角点扩展为8节点 (Q8 Serendipity)
% 输入：[x1, x2, x3, x4] - 角点顺序通常为逆时针
% 输出：[x1, x2, x3, x4, x5, x6, x7, x8] - 角点 + 边中点

x1 = xin(1); y1 = yin(1);  % 角点1
x2 = xin(2); y2 = yin(2);  % 角点2
x3 = xin(3); y3 = yin(3);  % 角点3
x4 = xin(4); y4 = yin(4);  % 角点4

% 计算边中点
x5 = (x1 + x2) / 2;  y5 = (y1 + y2) / 2;  % 边1-2中点
x6 = (x2 + x3) / 2;  y6 = (y2 + y3) / 2;  % 边2-3中点
x7 = (x3 + x4) / 2;  y7 = (y3 + y4) / 2;  % 边3-4中点
x8 = (x4 + x1) / 2;  y8 = (y4 + y1) / 2;  % 边4-1中点

% 8节点顺序：4角点 + 4边中点
xout = [x1, x2, x3, x4, x5, x6, x7, x8];
yout = [y1, y2, y3, y4, y5, y6, y7, y8];
end

function [xout, yout] = expand_to_5nodes_down(xin, yin)
% 将4角点扩展为5节点 ("Down" collapsed)
% 参考 Main_halfspace.m 中的 typei1=11 定义

x1 = xin(1); y1 = yin(1);
x2 = xin(2); y2 = yin(2);
x3 = xin(3); y3 = yin(3);
x4 = xin(4); y4 = yin(4);

% 添加边2-3的中点
x5 = (x2 + x3) / 2;  y5 = (y2 + y3) / 2;

% 节点顺序：[1, 2, 3, 5, 4]
xout = [x1, x2, x3, x5, x4];
yout = [y1, y2, y3, y5, y4];
end

function [xout, yout] = expand_to_5nodes_right(xin, yin)
% 将4角点扩展为5节点 ("Right" collapsed)
% 参考 Main_halfspace.m 中的 typei2=12 定义

x1 = xin(1); y1 = yin(1);
x2 = xin(2); y2 = yin(2);
x3 = xin(3); y3 = yin(3);
x4 = xin(4); y4 = yin(4);

% 添加边4-1的中点
x5 = (x4 + x1) / 2;  y5 = (y4 + y1) / 2;

% 节点顺序：[1, 2, 3, 4, 5]
xout = [x1, x2, x3, x4, x5];
yout = [y1, y2, y3, y4, y5];
end
