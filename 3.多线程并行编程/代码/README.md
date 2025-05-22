# PCFG 密码猜测项目并行化实验说明

## 主要变更内容

### 1. 新增并行猜测实现

- **guessing_pthread.cpp**
  - 使用 C++11 标准线程库（std::thread）和线程池实现了基于多线程的并行口令猜测生成。
  - 主要对 `PriorityQueue::Generate` 和 `PopNextBatchParallel` 进行了并行化处理。
  - 支持任务分块、线程池调度，提升大规模猜测时的吞吐量。

- **guessing_openmp.cpp**
  - 使用 OpenMP 实现了基于多核的并行口令猜测生成。
  - 主要对 `PriorityQueue::Generate` 进行了 `#pragma omp parallel for` 并行加速。
  - 自动负载均衡，适合多核 CPU 环境下的高效猜测生成。

### 2. PCFG.h 接口扩展

- 在 `PriorityQueue` 类中新增了批量并行接口声明：
  ```cpp
  void PopNextBatchParallel(int batch_pt_num);
  ```
  - 支持批量并行出队和猜测生成，便于在主程序中灵活选择串行或并行批量处理。

### 3. 相关主程序与编译说明

- `main.cpp`、`correctness_guess.cpp` 等主程序可通过切换链接不同的 guessing 源文件（如 guessing.cpp、guessing_pthread.cpp、guessing_openmp.cpp）来选择不同的并行实现。
- 编译示例：
  ```sh
  g++ main.cpp train.cpp guessing.cpp md5.cpp -o main
  g++ main.cpp train.cpp guessing_pthread.cpp md5.cpp -o main -O2
  g++ main.cpp train.cpp guessing_openmp.cpp md5.cpp -o main -O2 -fopenmp
  ```

## 使用建议

- 推荐在多核机器上使用 guessing_openmp 或 guessing_pthread 版本以获得更高的猜测生成速度。
- 若需批量并行处理，可调用 `PopNextBatchParallel` 接口。
- 相关并行代码已考虑负载均衡和线程安全，适合大规模实验和性能测试。

## 其他说明

- 训练过程（train.cpp）仍为串行实现，保证模型统计的准确性和可复现性。
- 并行猜测部分已通过日志和实验验证负载均衡和正确性。
- 如需进一步优化训练或哈希部分，可参考 README 顶部的优化建议。

