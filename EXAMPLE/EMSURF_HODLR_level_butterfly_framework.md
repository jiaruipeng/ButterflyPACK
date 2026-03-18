# 从 `EMSURF_Driver_sp.f90` 出发的 HODLR (`level_butterfly > 0`) 算法框架

> 适用前提：`EMSURF_Driver_sp.f90` + `EMSURF_Module_sp.f90`，单复数配置 `#define DAT 2`。

## 1. 主调用链（Driver 视角）

在 `EXAMPLE/EMSURF_Driver_sp.f90` 中，HODLR 路径由以下顺序触发：

1. 设置 `option%format = HODLR`。
2. `c_BPACK_construction_Init(...)`：初始化树结构/分块元数据。
3. `c_BPACK_construction_Element(...)`：执行压缩构造。
4. `c_BPACK_Factorization(...)`：执行分解（预条件/求解前处理）。
5. `EM_solve_SURF(...)`：对激励右端项求解。

可以把它理解为：

```text
几何+核函数注册
    -> HODLR结构化
    -> HODLR构造(压缩)
    -> HODLR分解
    -> 迭代/直接求解
```

---

## 2. `level_butterfly` 的判定与含义

对 HODLR 的某个 off-diagonal block：

- 若达到低秩阈值（例如 `Maxlevel - level < LRlevel`），则 `level_butterfly = 0`，走纯 LR。
- 否则 `level_butterfly > 0`，走 butterfly 压缩。

在 HODLR 代码中，常见赋值形式为：

```text
level_butterfly = Maxlevel - block_level
```

因此：

- `=0`：矩阵块直接使用 `U*V^T` 低秩因子。
- `>0`：矩阵块需要多层 butterfly skeleton/kernel 结构。

---

## 3. 结构化阶段（Structuring）框架

`BPACK_structuring -> HODLR_structuring` 里完成：

1. 建立每层 `levels(level_c)` 的 forward/inverse block 容器。
2. 为每个对角块 (`BP_inverse`) 和非对角块 (`BP`) 填充：
   - `row_group/col_group`
   - `headm/headn`, `M/N`
   - `pgno` 与并行局部分块索引
3. 对非对角块计算 `block_f%level_butterfly`：
   - `0`：低秩块
   - `>0`：butterfly 块
4. 当 `level_butterfly > 0` 时，继续构建 `LL` 层次边界映射（boundary map），为后续 butterfly 随机构造准备通信/索引关系。

这一阶段的本质：**先定义“块长什么样、在哪个进程、采用 LR 还是 BF”**。

---

## 4. 构造阶段（Construction）框架：`level_butterfly > 0`

`BPACK_construction_Element -> HODLR_construction` 中，对每个 HODLR 非对角块：

1. 初始化随机采样参数（`rank0`, `rankrate`, 容差）。
2. 若 `level_butterfly == 0`：走低秩构造分支（常规 LR/ACA 路线）。
3. 若 `level_butterfly > 0`：调用 `BF_randomized(...)`。

`BF_randomized` 在逻辑上完成：

- 通过用户核函数 matvec（或元素采样）做随机投影；
- 逐层构建 butterfly 的 skeleton/kernel；
- 估计误差并按 `rankrate` 扩展秩，直到达到容差；
- 形成可用于 MVP/分解的 BF 表示。

可抽象为：

```text
for 每个 off-diagonal block:
    if level_butterfly == 0:
        LR_compress(block)
    else:
        BF_randomized(block, level_butterfly, tol, rank0, rankrate)
```

---

## 5. 分解阶段（Factorization）框架：`level_butterfly > 0`

`BPACK_Factorization -> HODLR_factorization` 中，核心是层次 Schur 补更新与块运算。

当涉及 off-diagonal 更新（Add / Multiply / XLM / XUM）：

- `level_butterfly == 0`：优先走确定性 LR 快路径。
- `level_butterfly > 0`：再次调用 `BF_randomized(...)`，但这次是对“块运算结果”做压缩（不是最初 A_ij 的原始构造）。

典型调用语义包括：

- `'Multiply'`
- `'Add'`
- `'Add_Multiply'`
- `'XLM'`, `'XUM'`

所以分解期的关键点是：

> **BF 不仅用于初始构造，也用于分解过程中新产生的中间块压缩。**

---

## 6. 一页式总览（建议作为阅读代码顺序）

```text
EMSURF_Driver_sp
  ├─ 设置 option%format = HODLR
  ├─ construction_Init
  │    └─ BPACK_structuring
  │         └─ HODLR_structuring
  │              └─ 计算每个块 level_butterfly
  ├─ construction_Element
  │    └─ HODLR_construction
  │         └─ if level_butterfly>0: BF_randomized(原始块压缩)
  ├─ Factorization
  │    └─ HODLR_factorization
  │         └─ if level_butterfly>0: BF_randomized(更新块压缩)
  └─ EM_solve_SURF
```

---

## 7. 调参与排障建议（围绕 `level_butterfly>0`）

1. **确认确实进入 BF 分支**：检查 `LRlevel` 是否过大（过大可能让许多块退化成 `level_butterfly=0`）。
2. **先稳后快**：初次调试可降低 `rankrate`、收紧容差，优先看误差是否单调下降。
3. **观察层级统计**：对每层 rank/time 做对比，若某层异常膨胀，通常是该层几何分割与核可分性较弱。
4. **并行通信热点**：`level_butterfly` 越高，跨层通信与 boundary map 的开销越敏感。

