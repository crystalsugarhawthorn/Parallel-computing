#include "PCFG.h"
#include <cstdint>
#include <mpi.h>
#include <omp.h>
#include <cstring>
#include <numeric>
using namespace std;

// 优化: 缓存查找结果，避免在同一个循环内重复调用 Find* 函数
void PriorityQueue::CalProb(PT &pt)
{
    pt.prob = pt.preterm_prob; // 设置初始概率
    int index = 0;

    for (int idx : pt.curr_indices)
    {
        const auto& seg_content = pt.content[index];
        if (seg_content.type == 1)
        {
            // 优化: 只调用一次 FindLetter
            int seg_model_idx = m.FindLetter(seg_content);
            const auto& letter_seg = m.letters[seg_model_idx];
            if (letter_seg.total_freq > 0) { // 避免除以零
                pt.prob *= static_cast<float>(letter_seg.ordered_freqs[idx]) / letter_seg.total_freq;
            }
        }
        else if (seg_content.type == 2)
        {
            // 优化: 只调用一次 FindDigit
            int seg_model_idx = m.FindDigit(seg_content);
            const auto& digit_seg = m.digits[seg_model_idx];
            if (digit_seg.total_freq > 0) { // 避免除以零
                pt.prob *= static_cast<float>(digit_seg.ordered_freqs[idx]) / digit_seg.total_freq;
            }
        }
        else if (seg_content.type == 3)
        {
            // 优化: 只调用一次 FindSymbol
            int seg_model_idx = m.FindSymbol(seg_content);
            const auto& symbol_seg = m.symbols[seg_model_idx];
            if (symbol_seg.total_freq > 0) { // 避免除以零
                pt.prob *= static_cast<float>(symbol_seg.ordered_freqs[idx]) / symbol_seg.total_freq;
            }
        }
        index++;
    }
}

void PriorityQueue::init()
{
    // 优化: 使用 const 引用避免在循环中对 m.ordered_pts 中的每个 PT 对象进行不必要的拷贝
    for (const PT& pt_template : m.ordered_pts)
    {
        PT pt = pt_template; // 创建一个可修改的副本
        pt.max_indices.clear(); // 确保 max_indices 是空的

        for (const segment& seg : pt.content)
        {
            if (seg.type == 1)
            {
                pt.max_indices.emplace_back(m.letters[m.FindLetter(seg)].ordered_values.size());
            }
            else if (seg.type == 2)
            {
                pt.max_indices.emplace_back(m.digits[m.FindDigit(seg)].ordered_values.size());
            }
            else if (seg.type == 3)
            {
                pt.max_indices.emplace_back(m.symbols[m.FindSymbol(seg)].ordered_values.size());
            }
        }

        if (m.total_preterm > 0) {
            pt.preterm_prob = static_cast<float>(m.preterm_freq[m.FindPT(pt)]) / m.total_preterm;
        } else {
            pt.preterm_prob = 0.0f;
        }

        CalProb(pt);
        // 修正: 使用 .push() 方法向 std::priority_queue 添加元素
        priority.push(pt);
    }
}

// MPI 版本的 PopNext 也进行同样的修正
void PriorityQueue::PopNext_MPI() {
    if (priority.empty()) return;

    // 修正: 使用 .top() 获取最高优先级的元素
    const PT top_pt = priority.top();
    
    // 修正: 使用 .pop() 移除顶部元素
    priority.pop();

    // 注: Generate 函数内部处理 MPI 通信，这里不用动
    // Generate(top_pt); // 在MPI场景下，Generate通常由root进程调用并广播，或者每个进程独立计算

    // 修正: 与 PopNext() 逻辑保持一致
    PT pt_to_expand = top_pt;
    vector<PT> new_pts = pt_to_expand.NewPTs();

    for (PT& new_pt : new_pts) {
        CalProb(new_pt);
        // 修正: 直接 push，优先队列会自动排序
        priority.push(std::move(new_pt));
    }
}

// PT::NewPTs 函数本身逻辑没有问题，无需修改
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
        // 优化: curr_indices.size() 是不变量，可以在循环外计算
        // 但 size() 通常是O(1)的，现代编译器优化后性能差异可忽略，保持原样可读性更好
        for (int i = pivot; i < curr_indices.size(); i += 1)
        {
            // 如果是最后一个 segment，它没有 value 可以实例化，是用来生成猜测的
            // 所以我们只对 size-1 之前的 segment 进行扩展
            if (i == content.size() - 1) continue;

            curr_indices[i] += 1;
            if (curr_indices[i] < max_indices[i])
            {
                pivot = i;
                res.emplace_back(*this);
            }
            curr_indices[i] -= 1; // 恢复状态，为下一个可能的扩展做准备
        }
        pivot = init_pivot; // 恢复 pivot
        return res;
    }
    return res;
}

// Generate 函数的并行逻辑看起来是合理的，主要优化集中在数据准备和通信上
// 当前的 MPI_Gatherv 实现已经是一种比较高效的批量通信方式，无需大改
// 在你的 guessing_mpi.cpp 或类似文件中修改此函数

// 修正: Generate 函数现在接收一个 vector 的引用，用于存放当前进程生成的本地猜测
// 它不再直接操作成员变量 q.guesses 或执行 MPI 通信
void PriorityQueue::Generate(const PT& pt, vector<string>& local_guesses) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    local_guesses.clear(); // 确保从一个干净的状态开始

    if (pt.content.size() == 1) {
        segment* a = nullptr;
        if (pt.content[0].type == 1) a = &m.letters[m.FindLetter(pt.content[0])];
        else if (pt.content[0].type == 2) a = &m.digits[m.FindDigit(pt.content[0])];
        else a = &m.symbols[m.FindSymbol(pt.content[0])];

        if (a && !pt.max_indices.empty()) {
            int total_vals = pt.max_indices[0];
            // 优化: 预分配内存
            local_guesses.reserve(total_vals / size + 1);
            for (int i = rank; i < total_vals; i += size) {
                local_guesses.emplace_back(a->ordered_values[i]);
            }
        }
    } else {
        string prefix;
        // ... (生成 prefix 的逻辑保持不变) ...
        int seg_idx = 0;
        for (int idx : pt.curr_indices) {
            const auto& seg_content = pt.content[seg_idx];
            if (seg_content.type == 1) prefix += m.letters[m.FindLetter(seg_content)].ordered_values[idx];
            else if (seg_content.type == 2) prefix += m.digits[m.FindDigit(seg_content)].ordered_values[idx];
            else prefix += m.symbols[m.FindSymbol(seg_content)].ordered_values[idx];
            seg_idx++;
        }
        
        segment* a = nullptr;
        int last = pt.content.size() - 1;
        if (pt.content[last].type == 1) a = &m.letters[m.FindLetter(pt.content[last])];
        else if (pt.content[last].type == 2) a = &m.digits[m.FindDigit(pt.content[last])];
        else a = &m.symbols[m.FindSymbol(pt.content[last])];

        if (a && !pt.max_indices.empty()) {
            int total_vals = pt.max_indices.back();
            // 优化: 预分配内存
            local_guesses.reserve(total_vals / size + 1);
            for (int i = rank; i < total_vals; i += size) {
                local_guesses.emplace_back(prefix + a->ordered_values[i]);
            }
        }
    }
}