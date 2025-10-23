# 有限块体法 - 线弹性求解器 (MATLAB)

基于强形式的有限块体法（Finite Block Method）求解线弹性力学问题的 MATLAB 实现。

## 描述

本项目实现了一个强形式配点法求解器，用于求解平面应变/应力和轴对称条件下的线弹性问题。该方法采用：

- **等参映射**：支持多种单元类型（Q8 Serendipity、四节点无限单元、退化单元）
- **张量积 Lagrange 基函数**：高精度插值和导数计算
- **强形式配点法**：直接在配点上满足 Navier 方程和边界条件
- **多块域离散**：支持复杂几何和界面连续性条件
- **JSON 输入文件**：灵活的算例配置
- **轴对称问题**：支持圆柱、旋转壳等轴对称结构（含 r→0 数值稳定性处理）

示例应用场景：
- 悬臂梁变形分析
- 含孔平板应力集中
- 多材料界面问题
- 半无限域问题
- 厚壁圆柱/压力容器（轴对称）

## 目录

- [安装](#安装)
- [使用方法](#使用方法)
- [功能特性](#功能特性)
- [配置说明](#配置说明)
- [许可证](#许可证)

## 安装

### 前置要求

```matlab
% MATLAB 版本要求
MATLAB R2019b 或更高版本

% 必需工具箱
- 无（仅使用 MATLAB 核心功能）
```

### 安装步骤

1. 克隆仓库

```bash
git clone https://github.com/hhwon/FiniteBlockMethod-LED-MATLAB.git
cd FiniteBlockMethod-LED-MATLAB
```

2. 在 MATLAB 中打开项目文件夹

```matlab
cd('path/to/FiniteBlockMethod-LED-MATLAB')
```

3. 验证安装

```matlab
% 运行测试算例
main
```

## 使用方法

### 快速开始

```matlab
% 1. 打开 MATLAB 并进入项目目录
cd('path/to/FiniteBlockMethod-LED-MATLAB')

% 2. 运行主程序（使用默认算例）
main

% 3. 查看结果
% 程序会自动生成可视化图形
```

### 基本工作流程

#### 步骤 1: 准备输入文件

创建 JSON 格式的算例文件，定义：
- 几何（单元节点坐标）
- 材料参数（杨氏模量、泊松比）
- 边界条件（位移约束、牵引力）
- 离散化参数（每个方向的配点数）

示例见 [`case1_cantilever_oneblock.json`](case1_cantilever_oneblock.json)

#### 步骤 2: 修改主程序

在 [`main.m`](main.m) 中指定输入文件：

```matlab
CASE = load_json_case('case1_cantilever_oneblock.json');
```

#### 步骤 3: 运行求解

```matlab
main
```

#### 步骤 4: 查看结果

程序会生成以下可视化图形：
- 单元节点分布
- 配点分布
- 算子稀疏结构
- 位移场云图（u 方向、v 方向、位移幅值、变形构型）

### 高级用法

#### 自定义算例

```matlab
% 创建自定义算例结构
CASE = struct();
CASE.nt = 21;  % 每方向配点数
CASE.blocks = struct('xcoor', {[0 1 1 0]}, 'ycoor', {[0 0 1 1]}, 'type', {0});
CASE.E = {1.0};
CASE.v = {0.3};
CASE.plane_strain = true;
CASE.bcs = [];  % 边界条件数组

% 运行求解流程
[xcoor, ycoor, type] = expand_nodes(CASE.blocks);
% ... 后续步骤见 main.m
```

#### 单独调用关键函数

```matlab
% 1. 生成配点
[xx, yy, xt, yt, span] = generate_collocation_points(xcoor, ycoor, type, 21);

% 2. 计算形函数
[Phi, Dphix, Dphiy, Dphixx, Dphixy, Dphiyy, delta11, delta12, delta21, delta22] = ...
    compute_shape_functions(xcoor, ycoor, type, xt, yt, span);

% 3. 装配刚度矩阵
[K, f] = assemble_global_stiffness(CASE, Aphi, ADphix, ADphiy, ...
                                   ADphixx, ADphixy, ADphiyy, ...
                                   xx, yy, xt, yt, nblock, nnodes, ...
                                   delta11, delta12, delta21, delta22);

% 4. 求解
u = K \ f;

% 5. 可视化
visualize_solution(u, xx, yy, nblock, nnodes);
```

## 功能特性

### 核心功能

- ✨ **多种单元类型**
  - `type=0`: Q8 Serendipity（8 节点）
  - `type=2`: 4 节点无限单元
  - `type=11/12`: 单向无限单元（5 节点）
  
- 🚀 **高效数值算法**
  - 重心插值公式（Lagrange 基函数）
  - 张量积构造二维基函数
  - 稀疏矩阵存储和求解
  
- 💡 **灵活的边界条件**
  - Dirichlet 边界（指定位移）
  - Neumann 边界（指定牵引力）
  - 自由边界（零应力）
  
- 🔧 **界面连续性**
  - 位移连续条件
  - 应力平衡条件
  - 支持多块域离散
  
- 📦 **完整的可视化**
  - 节点和配点分布
  - 算子稀疏结构
  - 位移场云图（contourf）
  - 变形构型（带放大系数）

### 支持的算例类型

1. **悬臂梁**（单块）
   - 文件: [`case1_cantilever_oneblock.json`](case1_cantilever_oneblock.json)
   - 左端固定，右端施加向下载荷

2. **悬臂梁**（双块）
   - 文件: [`case2_cantilever_twoblock.json`](case2_cantilever_twoblock.json)
   - 验证界面连续性条件

3. **含孔平板**（双块）
   - 文件: [`case3_platewithhole_twoblock.json`](case3_platewithhole_twoblock.json)
   - 曲边单元、应力集中

4. **厚壁圆柱**（轴对称）
   - 文件: [`case4_axisym_cylinder.json`](case4_axisym_cylinder.json)
   - 轴对称问题，内压载荷，可与 Lamé 解析解对比

### 轴对称模型特性

- ✅ **轴对称实现**
  - 支持 r-z 坐标系（径向-轴向）
  - 自动处理周向应变 εθθ = u/r
  - 修改后的 Navier 方程包含 1/r 项
  
- 🔬 **数值稳定性保证**
  - r→0 时使用 L'Hospital 法则：φ/r → ∂φ/∂r
  - 自动检测对称轴（r < 1e-10）
  - 避免除零错误
  
- 📊 **解析解验证**
  - 厚壁圆柱 Lamé 解
  - 内压/外压载荷
  - 应力和位移对比

## 配置说明

### JSON 输入文件格式

```json
{
  "geometry": {
    "blocks": [
      {
        "id": 0,
        "vertices": [[x1,y1], [x2,y2], [x3,y3], [x4,y4]],
        "type": 0
      }
    ]
  },
  "discretization": {
    "node": 21
  },
  "material": {
    "model": "linear_elastic",
    "plane_strain": true,
    "axisymmetric": false,
    "defaults": {
      "E": 1.0,
      "nu": 0.3
    }
  },
  "boundaryConditions": [
    {
      "blk": 0,
      "face": "left",
      "bc": "dirichlet_ux",
      "value": 0.0
    }
  ]
}
```

### 关键参数说明

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `discretization.node` | 每方向配点数（推荐奇数） | 21 |
| `material.E` | 杨氏模量 | 1.0 |
| `material.nu` | 泊松比 | 0.3 |
| `material.plane_strain` | 平面应变（`true`）或平面应力（`false`） | `true` |
| `material.axisymmetric` | 轴对称问题（`true`）或平面问题（`false`） | `false` |
| `geometry.blocks[].type` | 单元类型（0/2/11/12） | 0 |

**重要说明**：
- `plane_strain` 和 `axisymmetric` 均放在 `material` 字段中，因为它们影响本构关系（应力-应变关系）
- 轴对称问题中，坐标解释为 (r, z)，位移为 (u_r, u_z)，其中 r > 0
- 轴对称实现包含 r→0 时的数值稳定性处理：φ/r → ∂φ/∂r

### 边界条件类型

| 边界条件 | 说明 | 平面问题 | 轴对称问题 |
|----------|------|----------|-----------|
| `dirichlet_ux` | 指定 x/r 方向位移 | x 方向 | r 方向（径向） |
| `dirichlet_uy` | 指定 y/z 方向位移 | y 方向 | z 方向（轴向） |
| `traction_x` | 指定 x/r 方向牵引力 | x 方向 | r 方向（径向） |
| `traction_y` | 指定 y/z 方向牵引力 | y 方向 | z 方向（轴向） |

**注意**：在轴对称问题中，边界条件的含义自动转换为径向和轴向。

### 界面条件

```json
"interfaces": [
  {
    "master_blk": 0,
    "master_face": "right",
    "slave_blk": 1,
    "slave_face": "left"
  }
]
```

- **主侧**：施加应力平衡条件
- **从侧**：施加位移连续条件

## 项目结构

```
FiniteBlockMethod-LED-MATLAB/
├── main.m                              # 主程序入口
├── load_json_case.m                    # JSON 输入文件解析
├── expand_nodes.m                      # 单元节点扩展
├── generate_collocation_points.m       # 配点生成
├── compute_shape_functions.m           # 形函数计算
├── assemble_global_operators.m         # 全局算子装配
├── assemble_global_stiffness.m         # 刚度矩阵装配
├── assemble_strong_form_operator.m     # 强形式算子
├── assemble_traction_bc_operator.m     # 牵引力边界算子
├── compute_boundary_normal.m           # 边界法向量计算
├── mappings.m                          # 等参映射
├── larg.m                              # 1D Lagrange 基函数
├── diffLarg.m                          # 2D 导数算子
├── visualize_*.m                       # 可视化函数
├── print_node_summary.m                # 节点摘要打印
├── verify_axisymmetric_solution.m      # 轴对称解析解验证
├── case*.json                          # 算例输入文件
└── README.md                           # 本文档
```

## 核心算法说明

### 1. 等参映射

使用 [`mappings.m`](mappings.m) 实现参考空间到物理空间的映射：

$$
\mathbf{x}(\xi, \eta) = \sum_{i=1}^{n} N_i(\xi, \eta) \mathbf{x}_i
$$

### 2. 形函数构造

张量积 Lagrange 基函数：

$$
\phi(\xi, \eta) = L_x(\xi) \otimes L_y(\eta)
$$

### 3. 强形式 Navier 方程

#### 平面应变/应力

$$
(\lambda + 2\mu) \frac{\partial^2 u}{\partial x^2} + \mu \frac{\partial^2 u}{\partial y^2} + (\lambda + \mu) \frac{\partial^2 v}{\partial x \partial y} = 0
$$

$$
(\lambda + 2\mu) \frac{\partial^2 v}{\partial y^2} + \mu \frac{\partial^2 v}{\partial x^2} + (\lambda + \mu) \frac{\partial^2 u}{\partial x \partial y} = 0
$$

#### 轴对称（r-z 坐标）

$$
(\lambda + 2\mu) \frac{\partial^2 u}{\partial r^2} + \mu \frac{\partial^2 u}{\partial z^2} + (\lambda + \mu) \frac{\partial^2 w}{\partial r \partial z} + \frac{\lambda + 2\mu}{r} \frac{\partial u}{\partial r} - \frac{\lambda + 2\mu}{r^2} u = 0
$$

$$
(\lambda + 2\mu) \frac{\partial^2 w}{\partial z^2} + \mu \frac{\partial^2 w}{\partial r^2} + (\lambda + \mu) \frac{\partial^2 u}{\partial r \partial z} + \frac{\mu}{r} \frac{\partial w}{\partial r} = 0
$$

**关键差异**：轴对称方程包含额外的 $1/r$ 项，来自周向应变 $\varepsilon_{\theta\theta} = u/r$

### 4. 边界条件

- **Dirichlet**: $u = \bar{u}$ 或 $v = \bar{v}$（$w = \bar{w}$ 在轴对称中）
- **Neumann**: $\boldsymbol{\sigma} \cdot \mathbf{n} = \mathbf{t}$

### 5. 数值稳定性处理

在轴对称问题中，当 $r \to 0$ 时使用 L'Hôpital 法则：

$$
\lim_{r \to 0} \frac{\phi}{r} = \frac{\partial \phi}{\partial r}, \quad \lim_{r \to 0} \frac{u}{r} = \frac{\partial u}{\partial r}
$$

实现中使用阈值 `r_tol = 1e-5` 自动切换计算公式。


## 故障排除

### 常见问题

1. **刚度矩阵奇异**
   - 检查边界条件是否充分约束刚体运动
   - 确保至少有一个 Dirichlet 边界条件

2. **收敛性差**
   - 增加配点数（`discretization.node`）
   - 检查块体质量（避免过度扭曲）

3. **界面不连续**
   - 确保界面节点对齐
   - 检查主从面定义是否正确

4. **轴对称问题错误**
   - 确保所有 r 坐标 > 0（不包含对称轴本身，除非使用特殊处理）
   - 检查 `material.axisymmetric` 是否设为 `true`

## 许可证

本项目采用 [MIT](LICENSE) 许可证。

---

## 作者

- 作者：W. Huang
- 项目地址：https://github.com/hhwon/FiniteBlockMethod-LED-MATLAB


## 致谢

感谢所有为这个项目做出贡献的人。

---

## 快速链接

### 核心文件
- [主程序](main.m)
- [JSON 输入解析](load_json_case.m)
- [强形式算子](assemble_strong_form_operator.m)
- [边界条件算子](assemble_traction_bc_operator.m)

### 算例文件
- [算例 1：悬臂梁（单块）](case1_cantilever_oneblock.json)
- [算例 2：悬臂梁（双块）](case2_cantilever_twoblock.json)
- [算例 3：含孔平板](case3_platewithhole_twoblock.json)
- [算例 4：厚壁圆柱（轴对称）](case4_axisym_cylinder.json)

### 验证和可视化
- [轴对称解析解验证](verify_axisymmetric_solution.m)
- [可视化函数](visualize_solution.m)
