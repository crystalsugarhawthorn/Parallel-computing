# PCFG 密码猜测项目并行化实验说明

## 主要变更内容

### 1. 新增并行猜测实现

- **guessing_pthread.cpp**
  - 使用 C++11 标准线程库（std::thread）和线程池实现了基于多线程的并行口令猜测生成。
  - 主要对 `PriorityQueue::Generate` 和 `PopNextBatchParallel` 进行了并行化处理。
  - 支持任务分块、线程池调度，提升大规模猜测时的吞吐量。
  - 适合多线程环境，能够充分利用 CPU 资源。

- **guessing_openmp.cpp**
  - 使用 OpenMP 实现了基于多核的并行口令猜测生成。
  - 主要对 `PriorityQueue::Generate` 和 `PopNextBatchParallel` 进行了 `#pragma omp parallel for` 并行加速。
  - 自动负载均衡，适合多核 CPU 环境下的高效猜测生成。
  - 代码简洁，易于维护，适合快速实验和性能测试。

### 2. PCFG.h 接口扩展

- 在 `PriorityQueue` 类中新增了批量并行接口声明：
  ```cpp
  void PopNextBatchParallel(int batch_pt_num);
  ```
  - 支持批量并行出队和猜测生成，便于在主程序中灵活选择串行或并行批量处理。
  - 通过批量处理减少线程调度开销，进一步提升性能。

- 引入 `PTComparator` 比较器，结合 `std::priority_queue` 实现了基于堆的候选 PT 管理：
  - 自动维护候选 PT 的概率降序排序。
  - 替代了原有的手动排序逻辑，简化代码并提升效率。

### 3. 相关主程序与编译说明

- `main.cpp`、`correctness_guess.cpp` 等主程序可通过切换链接不同的 guessing 源文件（如 guessing.cpp、guessing_pthread.cpp、guessing_openmp.cpp）来选择不同的并行实现。
- 编译示例：
  ```sh
  g++ main.cpp train.cpp guessing.cpp md5.cpp -o main
  g++ main.cpp train.cpp guessing_pthread.cpp md5.cpp -o main -O2
  g++ main.cpp train.cpp guessing_openmp.cpp md5.cpp -o main -O2 -fopenmp
  ```

### 4. 优化底层算法逻辑

- 使用 `std::priority_queue` 替代原有的 `vector`，实现基于堆的 beam-search 策略：
  - 每次从堆中取出概率最高的 PT 进行扩展。
  - 自动排序，无需手动维护候选 PT 的顺序。
  - 提升了代码可读性和扩展性，同时显著优化了性能。

- 修改了 `PriorityQueue::init`、`PopNext` 和 `PopNextBatchParallel` 方法：
  - `init` 方法将初始 PT 直接 `push` 入优先队列。
  - `PopNext` 方法从堆中取出最高概率的 PT，扩展后将新生成的 PT 重新 `push` 入堆中。
  - `PopNextBatchParallel` 方法支持批量从堆中取出 PT，并行扩展后统一 `push` 回堆中。

---

## 使用建议

- **多线程环境**：推荐使用 `guessing_pthread.cpp`，适合多线程 CPU，能够充分利用线程池的调度能力。
- **多核环境**：推荐使用 `guessing_openmp.cpp`，适合多核 CPU，代码简洁且性能优异。
- **批量处理**：若需批量并行处理，可调用 `PopNextBatchParallel` 接口，进一步提升吞吐量。
- **大规模实验**：建议在多核机器上运行，并使用 `-O2` 或更高优化级别编译代码。

---

## 其他说明

- 训练过程（train.cpp）仍为串行实现，保证模型统计的准确性和可复现性。
- 并行猜测部分已通过日志和实验验证负载均衡和正确性。
- 如需进一步优化训练或哈希部分，可参考 README 顶部的优化建议。

