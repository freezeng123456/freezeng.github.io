---
title: "GPU 如何工作：从 CUDA 生命周期到 SIMT"
description: 从 CPU–GPU 协作、CUDA 程序生命周期与执行映射理解 GPU 的吞吐、延迟和性能边界
lang: zh
translation: en/computer-science/knowledge-notes/how-gpus-work
tags:
  - GPU
  - CUDA
  - AI-Infra
  - Systems
---

> [!note] 适用范围
> 驱动命令格式、队列实现与调度细节会随 GPU 架构、驱动和 CUDA 版本变化。本页在这些位置采用公开文档支持的抽象模型，具体行为应以对应版本的源码与文档为准。

## 核心心智模型

GPU 不是由 CPU 逐条遥控的“大号算术单元”。应用通过 CUDA Runtime 和驱动准备内存、装载设备代码，并把工作异步提交给 GPU。GPU 前端解析命令，把 thread block 分配给可容纳它的 SM；SM 再以 warp 为调度单位执行指令。

![CPU 与 GPU 异构计算系统](assets/diagrams/gpu/zh/heterogeneous-system.svg)

这张图区分了三层并行：

1. **主机与设备并行**：异步提交后，CPU 可以继续执行；
2. **block 级并行**：不同 thread block 可以分布到不同 SM；
3. **warp 级并行**：一个 SM 在多个就绪 warp 之间切换，以隐藏等待延迟。

## 为什么 GPU 不一定让整个程序更快

原文用两个长度均为 $N=2^{30}$ 的 `float` 数组计算

$$
y_i \leftarrow x_i+y_i.
$$

GPU 版本在 NVIDIA T4 上使用 `blockSize = 256`，编译目标为 `sm_75`：

```cpp
int blockSize = 256;
int numBlocks = (N + blockSize - 1) / blockSize;
add<<<numBlocks, blockSize>>>(N, d_x, d_y);
```

原文记录的单次结果如下：

| 版本            | 核心计算时间 | 程序 `real` 时间 | 误差 |
| --------------- | -----------: | ---------------: | ---: |
| 单线程 CPU 循环 |      3740 ms |         21.418 s |    0 |
| T4 GPU kernel   |   48.6738 ms |         19.413 s |    0 |

按原文的舍入比较，GPU kernel 约快 75 倍，但整个程序只缩短约 2 秒。数组分配、主机端初始化、H2D/D2H 传输、CUDA 上下文与同步都不包含在纯 kernel 时间里，却会进入端到端时间。

> [!important] 这不是通用的 CPU–GPU 跑分
> CPU 版本是串行循环，GPU 结果来自一次特定 T4 运行。它说明的是**局部 kernel 加速不等于端到端加速**。严谨比较还需固定预热、内存策略、编译优化、并行 CPU 基线、传输范围和重复次数。

一个实用的分解是

$$
T_{\text{end-to-end}}
=
T_{\text{setup}}
+T_{\text{H2D}}
+T_{\text{kernel}}
+T_{\text{D2H}}
+T_{\text{sync}}
+T_{\text{host}}.
$$

只有当可并行计算足够大、数据复用足够多，或数据本来就在 GPU 上时，较小的 $T_{\text{kernel}}$ 才更容易转化为端到端收益。

## CUDA 程序从编译到执行

### 1. Fat binary

`nvcc` 把主机代码交给主机编译器，同时处理 `__global__` 等设备代码。可执行文件可以包含：

- **SASS/cubin**：面向某个具体 GPU 架构的机器码，例如 `sm_75`；
- **PTX**：虚拟指令集表示，可由驱动 JIT 编译到目标架构；
- **host object**：CPU 上运行的普通 C/C++ 代码。

若 fat binary 中存在兼容的 cubin，驱动可以直接装载；否则可以从 PTX JIT 编译，并通常缓存编译结果。精确的兼容选择受工具链和驱动版本约束。

### 2. 首次 CUDA 调用

进程第一次使用 CUDA 时，运行时和驱动需要建立设备上下文、地址空间与必要的管理状态。首次调用因此常比后续调用慢。性能测量应显式预热，避免把一次性初始化归入稳态 kernel 延迟。

### 3. Kernel launch

Kernel launch 会提交入口地址、Grid/Block 维度、动态共享内存大小和参数地址等信息。可以用 command queue 与 doorbell 建立一个足够准确的抽象：

![一次 CUDA Kernel Launch 抵达 SM 的过程](assets/diagrams/gpu/zh/kernel-launch.svg)

“Doorbell”表示 CPU 通过设备可见的寄存器通知 GPU 有新工作。具体命令包、队列位置和取数方式是驱动与硬件实现细节，不应依赖单一架构的内部描述编写应用逻辑。

`cudaMalloc` 传统上是同步、代价较高的分配路径；CUDA 11.2 起的 stream-ordered allocator 提供 `cudaMallocAsync` / `cudaFreeAsync`。工程上应优先复用内存池，避免在热路径频繁分配和释放。

## 编程模型如何映射到硬件

![CUDA 编程模型到 GPU 硬件执行模型的映射](assets/diagrams/gpu/zh/programming-model.svg)

| CUDA 抽象        | 硬件含义与约束                                                                   |
| ---------------- | -------------------------------------------------------------------------------- |
| Grid             | 一次 kernel launch 的全部 thread block                                           |
| Thread block     | 资源分配单位；放入某个 SM 后通常驻留到完成，block 内线程可共享内存并执行 barrier |
| Warp             | NVIDIA GPU 的基本调度组，通常包含 32 个线程                                      |
| Thread           | 独立的逻辑执行上下文，拥有索引和寄存器状态                                       |
| SM               | 容纳多个 block/warp，并含 warp scheduler、寄存器文件、共享内存和执行流水线       |
| CUDA/Tensor Core | 执行算术指令的流水线或专用单元，不等同于一个 CUDA thread                         |

一维数组常用

$$
i=\texttt{blockIdx.x}\times\texttt{blockDim.x}+\texttt{threadIdx.x}
$$

把软件线程映射到元素。相邻 warp lane 访问相邻地址时，内存请求更容易合并。二维 Grid/Block 则让矩阵与图像的坐标、tile 和数据局部性更自然。

## 内存层次与数据复用

GPU 性能不只由运算单元数量决定。程序还要决定数据在哪一级存储中停留：

| 层级                | 作用范围                     | 典型用途                                  |
| ------------------- | ---------------------------- | ----------------------------------------- |
| Registers           | 每线程                       | 标量、地址和局部中间量                    |
| Shared memory / L1  | 每 SM；共享内存由 block 使用 | tile 复用、block 内协作、降低全局内存流量 |
| L2                  | 全 GPU                       | 跨 SM 缓存与全局内存访问汇聚              |
| HBM / device memory | 全 GPU                       | 模型权重、激活、KV cache 和大数组         |
| Host memory         | CPU 侧                       | 通过 PCIe/NVLink-C2C 等链路与设备交换数据 |

关键不是背诵一个固定的延迟排序，而是减少高层级字节移动：合并全局内存访问、在 shared memory 或寄存器中复用 tile，并避免不必要的 H2D/D2H 往返。

## Occupancy 如何隐藏延迟

当一个 warp 等待内存或数据依赖时，warp scheduler 可以发射另一个就绪 warp。SM 上需要有足够多的 resident warps，但能驻留多少 block 同时受到寄存器、共享内存、线程数和架构上限约束。忽略资源分配粒度时，可用下面的简化上界建立直觉：

$$
B_{\mathrm{resident}}
\le
\min\left(
\left\lfloor\frac{R_{\mathrm{SM}}}{R_{\mathrm{thread}}T_{\mathrm{block}}}\right\rfloor,
\left\lfloor\frac{S_{\mathrm{SM}}}{S_{\mathrm{block}}}\right\rfloor,
\left\lfloor\frac{T_{\mathrm{SM}}}{T_{\mathrm{block}}}\right\rfloor,
B_{\mathrm{arch}}
\right).
$$

原文的寄存器示例中，一个 SM 有 65536 个寄存器，最多 2048 个线程。若每线程使用 64 个寄存器、每 block 256 个线程，则寄存器最多容纳 4 个 block，即 1024 个线程或 32 个 warp，对应 50% occupancy。若每线程增至 128 个寄存器，寄存器上限降为 2 个 block。

共享内存也会成为瓶颈：若每个 SM 有 96 KB、每 block 使用 32 KB，则最多同时驻留 3 个此类 block。

> [!warning] Occupancy 是约束，不是目标函数
> 100% occupancy 也可能全部 warp 都在等待内存。反过来，较低 occupancy 的计算密集 kernel 仍可能充分利用执行流水线。优化时应把 occupancy 与指令吞吐、内存带宽、缓存命中率和算术强度一起观察，可借助 Nsight Compute 与 roofline 分析定位瓶颈。

## SIMD、SIMT 与分支发散

| 模型 | 程序员看到的抽象                  | 硬件执行重点                                  |
| ---- | --------------------------------- | --------------------------------------------- |
| SIMD | 一条向量指令显式操作多个数据 lane | 程序或编译器负责构造向量与掩码                |
| SIMT | 编写单个 thread 的标量程序        | 硬件把线程组成 warp，并用活动掩码执行共同指令 |

SIMT 降低了显式向量化的编程负担，但没有消除 SIMD 式执行约束。如果同一 warp 的线程在 `if/else` 中选择不同路径，硬件通常需要分时执行各路径，并对不参与当前路径的 lane 关闭活动掩码：

![SIMT Warp Divergence 示意](assets/diagrams/gpu/zh/warp-divergence.svg)

Volta 及后续架构引入 Independent Thread Scheduling，保存更细粒度的线程执行状态并允许更灵活的重汇合与同步。它不意味着同一 warp 的 32 个线程可以在同一周期任意执行 32 条不同指令，warp divergence 的吞吐损失仍然存在。

同步原语的范围必须区分：

- `__syncthreads()` 是 block 级 barrier；所有应参与的未退出线程必须以一致方式到达；
- `__syncwarp(mask)` 只同步掩码指定的 warp lane，不一定是全部 32 个线程；
- 不同 block 之间不能用 `__syncthreads()` 同步，需要拆分 kernel、使用 cooperative groups，或借助其他全局协调机制。

## 性能诊断清单

1. **先分离时间**：初始化、内存分配、H2D/D2H、kernel、同步和 CPU 后处理分别计时。
2. **先预热再测量**：排除上下文创建、模块装载与 JIT 的一次性成本。
3. **检查访问模式**：相邻 lane 是否访问连续地址，tile 是否在 shared memory/寄存器中复用。
4. **检查资源限制**：registers/thread、shared memory/block、threads/block 分别允许多少 resident blocks。
5. **检查分支结构**：热点分支是否让同一 warp 长时间走不同路径。
6. **判断瓶颈类型**：用 roofline 与 profiler 区分 compute-bound、memory-bound 和 launch-bound。
7. **最后再调 block size**：`256` 是常见起点，不是普适最优值；应结合 kernel 资源用量与实测选择。

## 要点回顾

1. CPU 负责组织和提交工作，GPU 通过队列异步消费命令；同步点决定主机何时必须等待。
2. CUDA 的 Grid/Block/Thread 是软件抽象，SM/Warp/执行流水线是硬件执行层；block 与 warp 是连接两层的关键。
3. GPU kernel 可以比串行 CPU 循环快很多，但端到端收益还受初始化、传输、内存分配和主机工作支配。
4. Occupancy 提供可切换的 warp 来隐藏延迟，但高 occupancy 本身不保证高吞吐。
5. SIMT 隐藏了显式 SIMD 编程，无法消除分支发散；正确的同步范围与规整的访问模式同样重要。

## 主要参考资料

- [CUDA C++ Programming Guide](https://docs.nvidia.com/cuda/cuda-c-programming-guide/)
- [CUDA C++ Best Practices Guide](https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/)
- [Volta Tuning Guide: Independent Thread Scheduling](https://docs.nvidia.com/cuda/volta-tuning-guide/)
- [An Even Easier Introduction to CUDA](https://developer.nvidia.com/blog/even-easier-introduction-cuda/)
