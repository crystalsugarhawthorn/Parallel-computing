#include "PCFG.h"
#include <cstdint>
#include <cuda_runtime.h>
#include <vector>
#include <cstring>
#include <algorithm> // for std::make_heap, std::push_heap, std::pop_heap
using namespace std;

void PriorityQueue::CalProb(PT &pt)
{
    // 计算PriorityQueue里面一个PT的流程如下：
    // 1. 首先需要计算一个PT本身的概率。例如，L6S1的概率为0.15
    // 2. 需要注意的是，Queue里面的PT不是“纯粹的”PT，而是除了最后一个segment以外，全部被value实例化的PT
    // 3. 所以，对于L6S1而言，其在Queue里面的实际PT可能是123456S1，其中“123456”为L6的一个具体value。
    // 4. 这个时候就需要计算123456在L6中出现的概率了。假设123456在所有L6 segment中的概率为0.1，那么123456S1的概率就是0.1*0.15

    // 计算一个PT本身的概率。后续所有具体segment value的概率，直接累乘在这个初始概率值上
    pt.prob = pt.preterm_prob;

    // index: 标注当前segment在PT中的位置
    int index = 0;


    for (int idx : pt.curr_indices)
    {
        // pt.content[index].PrintSeg();
        if (pt.content[index].type == 1)
        {
            // 下面这行代码的意义：
            // pt.content[index]：目前需要计算概率的segment
            // m.FindLetter(seg): 找到一个letter segment在模型中的对应下标
            // m.letters[m.FindLetter(seg)]：一个letter segment在模型中对应的所有统计数据
            // m.letters[m.FindLetter(seg)].ordered_values：一个letter segment在模型中，所有value的总数目
            pt.prob *= m.letters[m.FindLetter(pt.content[index])].ordered_freqs[idx];
            pt.prob /= m.letters[m.FindLetter(pt.content[index])].total_freq;
            // cout << m.letters[m.FindLetter(pt.content[index])].ordered_freqs[idx] << endl;
            // cout << m.letters[m.FindLetter(pt.content[index])].total_freq << endl;
        }
        if (pt.content[index].type == 2)
        {
            pt.prob *= m.digits[m.FindDigit(pt.content[index])].ordered_freqs[idx];
            pt.prob /= m.digits[m.FindDigit(pt.content[index])].total_freq;
            // cout << m.digits[m.FindDigit(pt.content[index])].ordered_freqs[idx] << endl;
            // cout << m.digits[m.FindDigit(pt.content[index])].total_freq << endl;
        }
        if (pt.content[index].type == 3)
        {
            pt.prob *= m.symbols[m.FindSymbol(pt.content[index])].ordered_freqs[idx];
            pt.prob /= m.symbols[m.FindSymbol(pt.content[index])].total_freq;
            // cout << m.symbols[m.FindSymbol(pt.content[index])].ordered_freqs[idx] << endl;
            // cout << m.symbols[m.FindSymbol(pt.content[index])].total_freq << endl;
        }
        index += 1;
    }
    // cout << pt.prob << endl;
}

void PriorityQueue::init()
{
    // cout << m.ordered_pts.size() << endl;
    // 用所有可能的PT，按概率降序填满整个优先队列
    for (PT pt : m.ordered_pts)
    {
        for (segment seg : pt.content)
        {
            if (seg.type == 1)
            {
                // 下面这行代码的意义：
                // max_indices用来表示PT中各个segment的可能数目。例如，L6S1中，假设模型统计到了100个L6，那么L6对应的最大下标就是99
                // （但由于后面采用了"<"的比较关系，所以其实max_indices[0]=100）
                // m.FindLetter(seg): 找到一个letter segment在模型中的对应下标
                // m.letters[m.FindLetter(seg)]：一个letter segment在模型中对应的所有统计数据
                // m.letters[m.FindLetter(seg)].ordered_values：一个letter segment在模型中，所有value的总数目
                pt.max_indices.emplace_back(m.letters[m.FindLetter(seg)].ordered_values.size());
            }
            if (seg.type == 2)
            {
                pt.max_indices.emplace_back(m.digits[m.FindDigit(seg)].ordered_values.size());
            }
            if (seg.type == 3)
            {
                pt.max_indices.emplace_back(m.symbols[m.FindSymbol(seg)].ordered_values.size());
            }
        }
        pt.preterm_prob = float(m.preterm_freq[m.FindPT(pt)]) / m.total_preterm;
        // pt.PrintPT();
        // cout << " " << m.preterm_freq[m.FindPT(pt)] << " " << m.total_preterm << " " << pt.preterm_prob << endl;

        // 计算当前pt的概率
        CalProb(pt);
        priority.push_back(pt);
    }
    // 建堆
    std::make_heap(priority.begin(), priority.end(), PTComparator());
}

// PT::NewPTs 实现
vector<PT> PT::NewPTs()
{
    vector<PT> res;
    if (content.size() == 1)
    {
        return res;
    }
    else
    {
        int init_pivot = pivot;
        for (int i = pivot; i < curr_indices.size() - 1; i += 1)
        {
            curr_indices[i] += 1;
            if (curr_indices[i] < max_indices[i])
            {
                pivot = i;
                res.emplace_back(*this);
            }
            curr_indices[i] -= 1;
        }
        pivot = init_pivot;
        return res;
    }
}

void PriorityQueue::PopNext()
{
    // 堆顶元素
    if (priority.empty()) return;
    std::pop_heap(priority.begin(), priority.end(), PTComparator());
    PT pt = priority.back();
    priority.pop_back();

    // 生成猜测
    vector<string> local_guesses;
    Generate(pt, local_guesses);
    guesses.insert(guesses.end(), local_guesses.begin(), local_guesses.end());

    // 生成新PT并插入堆
    vector<PT> new_pts = pt.NewPTs();
    for (PT& new_pt : new_pts)
    {
        CalProb(new_pt);
        priority.push_back(new_pt);
        std::push_heap(priority.begin(), priority.end(), PTComparator());
    }
}

// 串行猜测生成
void PriorityQueue::Generate(const PT& pt, vector<string>& local_guesses)
{
    CalProb(const_cast<PT&>(pt));
    if (pt.content.size() == 1)
    {
        segment *a = nullptr;
        if (pt.content[0].type == 1) a = &m.letters[m.FindLetter(pt.content[0])];
        else if (pt.content[0].type == 2) a = &m.digits[m.FindDigit(pt.content[0])];
        else if (pt.content[0].type == 3) a = &m.symbols[m.FindSymbol(pt.content[0])];
        for (int i = 0; i < pt.max_indices[0]; i++)
            local_guesses.emplace_back(a->ordered_values[i]);
        total_guesses += pt.max_indices[0];
    }
    else
    {
        string guess;
        int seg_idx = 0;
        for (int idx : pt.curr_indices)
        {
            if (pt.content[seg_idx].type == 1)
                guess += m.letters[m.FindLetter(pt.content[seg_idx])].ordered_values[idx];
            else if (pt.content[seg_idx].type == 2)
                guess += m.digits[m.FindDigit(pt.content[seg_idx])].ordered_values[idx];
            else if (pt.content[seg_idx].type == 3)
                guess += m.symbols[m.FindSymbol(pt.content[seg_idx])].ordered_values[idx];
            seg_idx++;
            if (seg_idx == pt.content.size() - 1) break;
        }
        segment *a = nullptr;
        if (pt.content.back().type == 1) a = &m.letters[m.FindLetter(pt.content.back())];
        else if (pt.content.back().type == 2) a = &m.digits[m.FindDigit(pt.content.back())];
        else if (pt.content.back().type == 3) a = &m.symbols[m.FindSymbol(pt.content.back())];
        for (int i = 0; i < pt.max_indices.back(); i++)
            local_guesses.emplace_back(guess + a->ordered_values[i]);
        total_guesses += pt.max_indices.back();
    }
}

#ifdef __CUDACC__
#include <cuda_runtime.h>
#include <future>
#include <cub/cub.cuh>  // 用于高效的并行原语

// 配置参数
constexpr int WARP_SIZE = 32;
constexpr int MAX_GUESS_LEN = 512;  // 根据实际情况调整
constexpr int MIN_GPU_BATCH_SIZE = 10000;  // 只有足够大的批次才用GPU

// 内存对齐优化
inline int align_to(int n, int alignment) { 
    return ((n + alignment - 1) / alignment) * alignment; 
}

// 优化的二分查找工具
__device__ __forceinline__ int binary_search_pt(const int* cumulative_counts, int n, int idx) {
    int left = 0, right = n;
    while (left < right) {
        int mid = (left + right) >> 1;
        if (cumulative_counts[mid + 1] <= idx)
            left = mid + 1;
        else
            right = mid;
    }
    return left;
}

// 优化的GPU kernel：使用共享内存和二分查找
__global__ void GenerateGuessesOptimized(
    const char* __restrict__ prefixes,
    const char* __restrict__ values, 
    const int* __restrict__ prefix_lens,
    const int* __restrict__ value_lens,
    const int* __restrict__ cumulative_counts,
    int n_pts,
    int total_guesses,
    int max_guess_len,
    char* __restrict__ output)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= total_guesses) return;

    // 二分查找PT索引
    int pt_idx = binary_search_pt(cumulative_counts, n_pts, idx);
    int local_idx = idx - cumulative_counts[pt_idx];
    int prefix_len = prefix_lens[pt_idx];
    int value_len = value_lens[pt_idx];

    const char* prefix_ptr = prefixes + pt_idx * max_guess_len;
    const char* value_ptr = values + (cumulative_counts[pt_idx] + local_idx) * value_len;
    char* out_ptr = output + idx * max_guess_len;

    // 合并内存拷贝
    int copy_len = prefix_len + value_len;

    // 拷贝前缀
    for (int i = 0; i < prefix_len; ++i) out_ptr[i] = prefix_ptr[i];
    // 拷贝值
    for (int i = 0; i < value_len; ++i) out_ptr[prefix_len + i] = value_ptr[i];
    // 清零剩余部分
    for (int i = copy_len; i < max_guess_len; ++i) out_ptr[i] = 0;
}

class GPUMemoryPool {
private:
    struct Buffer {
        void* ptr;
        size_t size;
        bool in_use;
    };
    
    std::vector<Buffer> buffers;
    cudaStream_t stream;
    
public:
    GPUMemoryPool() {
        cudaStreamCreate(&stream);
    }
    
    ~GPUMemoryPool() {
        for (auto& buf : buffers) {
            cudaFree(buf.ptr);
        }
        cudaStreamDestroy(stream);
    }
    
    void* allocate(size_t size) {
        // 寻找合适的缓冲区或分配新的
        for (auto& buf : buffers) {
            if (!buf.in_use && buf.size >= size) {
                buf.in_use = true;
                return buf.ptr;
            }
        }
        
        // 分配新缓冲区
        void* ptr;
        cudaMalloc(&ptr, size);
        buffers.push_back({ptr, size, true});
        return ptr;
    }
    
    void deallocate(void* ptr) {
        for (auto& buf : buffers) {
            if (buf.ptr == ptr) {
                buf.in_use = false;
                break;
            }
        }
    }
    
    cudaStream_t get_stream() { return stream; }
};

// 单例内存池
GPUMemoryPool& get_memory_pool() {
    static GPUMemoryPool pool;
    return pool;
}

#endif

// 性能分析器
class GPUProfiler {
private:
    std::vector<std::pair<std::string, float>> timings;
    cudaEvent_t start, stop;
    
public:
    GPUProfiler() {
        cudaEventCreate(&start);
        cudaEventCreate(&stop);
    }
    
    ~GPUProfiler() {
        cudaEventDestroy(start);
        cudaEventDestroy(stop);
    }
    
    void start_timer(const std::string& name) {
        current_name = name;
        cudaEventRecord(start);
    }
    
    void end_timer() {
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        float ms;
        cudaEventElapsedTime(&ms, start, stop);
        timings.push_back({current_name, ms});
    }
    
    void print_stats() {
        std::cout << "GPU Performance Stats:\n";
        for (const auto& [name, time] : timings) {
            std::cout << name << ": " << time << " ms\n";
        }
    }
    
private:
    std::string current_name;
};

// 增加全局pinned memory池，避免频繁分配释放
struct GPUPinnedBuffer {
    char* prefixes = nullptr;
    char* values = nullptr;
    char* output = nullptr;
    size_t prefixes_size = 0;
    size_t values_size = 0;
    size_t output_size = 0;
    ~GPUPinnedBuffer() {
        if (prefixes) cudaFreeHost(prefixes);
        if (values) cudaFreeHost(values);
        if (output) cudaFreeHost(output);
    }
    void alloc(size_t pre_sz, size_t val_sz, size_t out_sz) {
        if (pre_sz > prefixes_size) {
            if (prefixes) cudaFreeHost(prefixes);
            cudaMallocHost(&prefixes, pre_sz);
            prefixes_size = pre_sz;
        }
        if (val_sz > values_size) {
            if (values) cudaFreeHost(values);
            cudaMallocHost(&values, val_sz);
            values_size = val_sz;
        }
        if (out_sz > output_size) {
            if (output) cudaFreeHost(output);
            cudaMallocHost(&output, out_sz);
            output_size = out_sz;
        }
    }
};
static GPUPinnedBuffer g_pinned_buf;

// 极致优化的GPU批量猜测生成
void PriorityQueue::Generate_GPU(const vector<PT>& pts, vector<string>& guesses_out) {
#ifdef __CUDACC__
    int n = pts.size();
    if (n == 0) return;

    // 1. 快速评估：如果总猜测数太少，直接用CPU
    int total_estimate = 0;
    for (const auto& pt : pts) {
        int count = (pt.content.size() == 1) ? pt.max_indices[0] : pt.max_indices.back();
        total_estimate += count;
        if (total_estimate > MIN_GPU_BATCH_SIZE) break;
    }
    if (total_estimate < MIN_GPU_BATCH_SIZE) {
        Generate_Serial_Optimized(pts, guesses_out);
        return;
    }

    GPUProfiler profiler;
    auto& pool = get_memory_pool();

    // 2. 预处理：计算所有必要的信息
    profiler.start_timer("Preprocessing");

    vector<string> all_prefixes, all_values;
    vector<int> prefix_lens, value_lens, value_counts, cumulative_counts(n + 1, 0);
    int max_guess_len = 0;

    for (int i = 0; i < n; ++i) {
        const auto& pt = pts[i];
        string prefix;
        if (pt.content.size() > 1) {
            int seg_idx = 0;
            for (int idx : pt.curr_indices) {
                if (pt.content[seg_idx].type == 1)
                    prefix += m.letters[m.FindLetter(pt.content[seg_idx])].ordered_values[idx];
                else if (pt.content[seg_idx].type == 2)
                    prefix += m.digits[m.FindDigit(pt.content[seg_idx])].ordered_values[idx];
                else
                    prefix += m.symbols[m.FindSymbol(pt.content[seg_idx])].ordered_values[idx];
                seg_idx++;
                if (seg_idx == pt.content.size() - 1) break;
            }
        }
        segment* last_seg = nullptr;
        if (pt.content.back().type == 1) 
            last_seg = &m.letters[m.FindLetter(pt.content.back())];
        else if (pt.content.back().type == 2) 
            last_seg = &m.digits[m.FindDigit(pt.content.back())];
        else 
            last_seg = &m.symbols[m.FindSymbol(pt.content.back())];
        int count = (pt.content.size() == 1) ? pt.max_indices[0] : pt.max_indices.back();
        int value_len = last_seg->ordered_values.empty() ? 0 : last_seg->ordered_values[0].size();
        all_prefixes.push_back(prefix);
        prefix_lens.push_back(prefix.size());
        value_lens.push_back(value_len);
        value_counts.push_back(count);
        cumulative_counts[i + 1] = cumulative_counts[i] + count;
        max_guess_len = max(max_guess_len, (int)prefix.size() + value_len);
        for (int j = 0; j < count; ++j) {
            all_values.push_back(last_seg->ordered_values[j]);
        }
    }
    int total_guesses = cumulative_counts[n];
    max_guess_len = align_to(max_guess_len, 4);

    profiler.end_timer();

    // 3. 高效的内存分配和数据传输
    profiler.start_timer("Memory Transfer");

    // 复用全局 pinned buffer
    g_pinned_buf.alloc(n * max_guess_len, total_guesses * max_guess_len, total_guesses * max_guess_len);
    char *h_prefixes_pinned = g_pinned_buf.prefixes;
    char *h_values_pinned = g_pinned_buf.values;
    char *h_output_pinned = g_pinned_buf.output;

    // 数据打包
    if (total_guesses > 100000) {
        #pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            memset(h_prefixes_pinned + i * max_guess_len, 0, max_guess_len);
            memcpy(h_prefixes_pinned + i * max_guess_len, all_prefixes[i].c_str(), all_prefixes[i].size());
        }
        #pragma omp parallel for
        for (int i = 0; i < total_guesses; ++i) {
            memset(h_values_pinned + i * max_guess_len, 0, max_guess_len);
            memcpy(h_values_pinned + i * max_guess_len, all_values[i].c_str(), all_values[i].size());
        }
    } else {
        for (int i = 0; i < n; ++i) {
            memset(h_prefixes_pinned + i * max_guess_len, 0, max_guess_len);
            memcpy(h_prefixes_pinned + i * max_guess_len, all_prefixes[i].c_str(), all_prefixes[i].size());
        }
        for (int i = 0; i < total_guesses; ++i) {
            memset(h_values_pinned + i * max_guess_len, 0, max_guess_len);
            memcpy(h_values_pinned + i * max_guess_len, all_values[i].c_str(), all_values[i].size());
        }
    }

    char* d_prefixes = (char*)pool.allocate(n * max_guess_len);
    char* d_values = (char*)pool.allocate(total_guesses * max_guess_len);
    char* d_output = (char*)pool.allocate(total_guesses * max_guess_len);
    int* d_prefix_lens = (int*)pool.allocate(n * sizeof(int));
    int* d_value_lens = (int*)pool.allocate(n * sizeof(int));
    int* d_cumulative_counts = (int*)pool.allocate((n + 1) * sizeof(int));

    cudaStream_t stream = pool.get_stream();
    cudaMemcpyAsync(d_prefixes, h_prefixes_pinned, n * max_guess_len, cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(d_values, h_values_pinned, total_guesses * max_guess_len, cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(d_prefix_lens, prefix_lens.data(), n * sizeof(int), cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(d_value_lens, value_lens.data(), n * sizeof(int), cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(d_cumulative_counts, cumulative_counts.data(), (n + 1) * sizeof(int), cudaMemcpyHostToDevice, stream);

    profiler.end_timer();

    // 4. 优化的kernel启动
    profiler.start_timer("Kernel Execution");

    int device_id;
    cudaGetDevice(&device_id);
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, device_id);
    int max_threads_per_block = prop.maxThreadsPerBlock;
    int block_size = min(1024, max_threads_per_block);
    int grid_size = (total_guesses + block_size - 1) / block_size;

    GenerateGuessesOptimized<<<grid_size, block_size, 0, stream>>>(
        d_prefixes, d_values, d_prefix_lens, d_value_lens, 
        d_cumulative_counts, n, total_guesses, max_guess_len, d_output);

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "Kernel launch failed: %s\n", cudaGetErrorString(err));
        return;
    }

    profiler.end_timer();

    // 5. 高效的结果回传
    profiler.start_timer("Result Transfer");

    cudaMemcpyAsync(h_output_pinned, d_output, total_guesses * max_guess_len, cudaMemcpyDeviceToHost, stream);
    cudaStreamSynchronize(stream);

    guesses_out.clear();
    guesses_out.reserve(total_guesses);

    // 并行构建结果字符串
    #pragma omp parallel for
    for (int i = 0; i < total_guesses; ++i) {
        int pt_idx = 0;
        // 二分查找pt_idx
        int l = 0, r = n;
        while (l < r) {
            int mid = (l + r) >> 1;
            if (cumulative_counts[mid + 1] <= i)
                l = mid + 1;
            else
                r = mid;
        }
        pt_idx = l;
        int actual_len = prefix_lens[pt_idx] + value_lens[pt_idx];
        guesses_out.emplace_back(h_output_pinned + i * max_guess_len, actual_len);
    }

    profiler.end_timer();

    // 6. 清理资源
    pool.deallocate(d_prefixes);
    pool.deallocate(d_values);
    pool.deallocate(d_output);
    pool.deallocate(d_prefix_lens);
    pool.deallocate(d_value_lens);
    pool.deallocate(d_cumulative_counts);

    // pinned buffer复用，不再cudaFreeHost

    // profiler.print_stats();  // 调试时启用

#else
    Generate_Serial_Optimized(pts, guesses_out);
#endif
}

// 优化的串行版本作为对比
void PriorityQueue::Generate_Serial_Optimized(const vector<PT>& pts, vector<string>& guesses_out) {
    guesses_out.clear();
    
    // 预分配内存
    int total_estimate = 0;
    for (const auto& pt : pts) {
        int count = (pt.content.size() == 1) ? pt.max_indices[0] : pt.max_indices.back();
        total_estimate += count;
    }
    guesses_out.reserve(total_estimate);
    
    // 优化的串行生成
    for (const auto& pt : pts) {
        string prefix;
        if (pt.content.size() > 1) {
            prefix.reserve(256);  // 预分配
            int seg_idx = 0;
            for (int idx : pt.curr_indices) {
                if (pt.content[seg_idx].type == 1)
                    prefix += m.letters[m.FindLetter(pt.content[seg_idx])].ordered_values[idx];
                else if (pt.content[seg_idx].type == 2)
                    prefix += m.digits[m.FindDigit(pt.content[seg_idx])].ordered_values[idx];
                else
                    prefix += m.symbols[m.FindSymbol(pt.content[seg_idx])].ordered_values[idx];
                seg_idx++;
                if (seg_idx == pt.content.size() - 1) break;
            }
        }
        
        segment* last_seg = nullptr;
        if (pt.content.back().type == 1) 
            last_seg = &m.letters[m.FindLetter(pt.content.back())];
        else if (pt.content.back().type == 2) 
            last_seg = &m.digits[m.FindDigit(pt.content.back())];
        else 
            last_seg = &m.symbols[m.FindSymbol(pt.content.back())];
        
        int count = (pt.content.size() == 1) ? pt.max_indices[0] : pt.max_indices.back();
        
        for (int i = 0; i < count; ++i) {
            guesses_out.emplace_back(prefix + last_seg->ordered_values[i]);
        }
    }
}

