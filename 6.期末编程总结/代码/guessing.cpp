<<<<<<< HEAD
#include "PCFG.h"
#include <queue>
#include <string_view>
#include <unordered_map>
#include <cmath>
#include <memory_resource>
#include <algorithm>
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
    size_t n = m.ordered_pts.size();
    vector<PT> pts(n);
    
    // 预计算 segment 索引，避免重复调用 Find* 函数
    vector<vector<int>> letter_indices(m.letters.size());
    vector<vector<int>> digit_indices(m.digits.size());
    vector<vector<int>> symbol_indices(m.symbols.size());
    
    #pragma omp parallel for
    for (size_t i = 0; i < m.letters.size(); ++i) {
        letter_indices[i].resize(m.letters[i].ordered_values.size());
        for (size_t j = 0; j < m.letters[i].ordered_values.size(); ++j) {
            letter_indices[i][j] = i; // 缓存索引以加快查找
        }
    }
    #pragma omp parallel for
    for (size_t i = 0; i < m.digits.size(); ++i) {
        digit_indices[i].resize(m.digits[i].ordered_values.size());
        for (size_t j = 0; j < m.digits[i].ordered_values.size(); ++j) {
            digit_indices[i][j] = i;
        }
    }
    #pragma omp parallel for
    for (size_t i = 0; i < m.symbols.size(); ++i) {
        symbol_indices[i].resize(m.symbols[i].ordered_values.size());
        for (size_t j = 0; j < m.symbols[i].ordered_values.size(); ++j) {
            symbol_indices[i][j] = i;
        }
    }

    // 并行初始化 PT 对象
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < n; ++i) {
        PT pt = m.ordered_pts[i];
        pt.max_indices.reserve(pt.content.size()); // 预分配 max_indices
        pt.max_indices.clear();
        
        for (const segment& seg : pt.content) {
            if (seg.type == 1) {
                int idx = m.FindLetter(seg);
                pt.max_indices.emplace_back(m.letters[idx].ordered_values.size());
            } else if (seg.type == 2) {
                int idx = m.FindDigit(seg);
                pt.max_indices.emplace_back(m.digits[idx].ordered_values.size());
            } else if (seg.type == 3) {
                int idx = m.FindSymbol(seg);
                pt.max_indices.emplace_back(m.symbols[idx].ordered_values.size());
            }
        }
        pt.preterm_prob = static_cast<float>(m.preterm_freq[m.FindPT(pt)]) / m.total_preterm;
        CalProb(pt);
        pts[i] = pt;
    }

    // 并行插入优先队列
    int n_threads = omp_get_max_threads();
    vector<vector<PT>> local_pts(n_threads);
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        local_pts[tid].reserve(n / n_threads + 1); // 每个线程预分配空间
        #pragma omp for schedule(static)
        for (size_t i = 0; i < n; ++i) {
            local_pts[tid].push_back(pts[i]);
        }
    }

    // 串行合并 local_pts 到优先队列
    for (const auto& thread_pts : local_pts) {
        for (const auto& pt : thread_pts) {
            priority.push(pt);
        }
    }
}

void PriorityQueue::PopNext()
{
    // 动态调整弹出数量，基于线程数和队列大小
    int n_threads = omp_get_max_threads();
    int pop_count = max(32, n_threads * 4); // 动态设置，至少32或线程数的4倍
    int actual_pop = min(pop_count, static_cast<int>(priority.size()));
    
    // 预分配 popped_pts 容器
    vector<PT> popped_pts;
    popped_pts.reserve(actual_pop);
    
    // 串行弹出 PT 对象（优先队列操作难以并行）
    for (int i = 0; i < actual_pop && !priority.empty(); ++i) {
        popped_pts.push_back(priority.top());
        priority.pop();
    }

    // 并行处理 Generate
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < popped_pts.size(); ++i) {
        //Generate(popped_pts[i]);
        Generate_Beam(popped_pts[i], 1024); // 使用新的 Generate_Beam 函数
    }

    // 并行生成新 PT 并收集到线程本地容器
    vector<vector<PT>> local_new_pts(n_threads);
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        local_new_pts[tid].reserve(popped_pts.size() / n_threads + 1); // 预分配本地容器
        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < popped_pts.size(); ++i) {
            vector<PT> new_pts = popped_pts[i].NewPTs();
            for (PT& pt : new_pts) {
                CalProb(pt);
                local_new_pts[tid].push_back(pt);
            }
        }
    }

    // 串行合并到优先队列
    for (const auto& thread_pts : local_new_pts) {
        for (const auto& pt : thread_pts) {
            priority.push(pt);
        }
    }
}

// 这个函数你就算看不懂，对并行算法的实现影响也不大
// 当然如果你想做一个基于多优先队列的并行算法，可能得稍微看一看了
vector<PT> PT::NewPTs()
{
    vector<PT> res;
    if (content.size() == 1)
        return res;

    int init_pivot = pivot;
    int n = curr_indices.size() - 1;
    res.reserve(n - init_pivot); // 预分配空间

    // 并行化收益有限，这里保持串行，保证顺序和低开销
    for (int i = init_pivot; i < n; ++i)
    {
        curr_indices[i] += 1;
        if (curr_indices[i] < max_indices[i])
        {
            pivot = i;
            res.emplace_back(*this); // 只在需要时拷贝
        }
        curr_indices[i] -= 1;
    }
    pivot = init_pivot;
    return res;
}


// 这个函数是PCFG并行化算法的主要载体
// 尽量看懂，然后进行并行实现
void PriorityQueue::Generate(PT pt)
{
    int max_total = 1000000; // 单个PT最大生成数量，建议根据实际需求调整
    // 计算PT的概率，这里主要是给PT的概率进行初始化
    CalProb(pt);

    // 对于只有一个segment的PT，直接遍历生成其中的所有value即可
    if (pt.content.size() == 1)
    {
        // 指向最后一个segment的指针，这个指针实际指向模型中的统计数据
        segment *a;
        // 在模型中定位到这个segment
        if (pt.content[0].type == 1)
        {
            a = &m.letters[m.FindLetter(pt.content[0])];
        }
        if (pt.content[0].type == 2)
        {
            a = &m.digits[m.FindDigit(pt.content[0])];
        }
        if (pt.content[0].type == 3)
        {
            a = &m.symbols[m.FindSymbol(pt.content[0])];
        }
        
        // Multi-thread TODO：
        // 这个for循环就是你需要进行并行化的主要部分了，特别是在多线程&GPU编程任务中
        // 可以看到，这个循环本质上就是把模型中一个segment的所有value，赋值到PT中，形成一系列新的猜测
        // 这个过程是可以高度并行化的
        int n = pt.max_indices[0];
        // n 就是当前PT的最大猜测数（即该segment的所有value数量）
        if (n > max_total)
        {
            n = max_total; // 限制单个PT的最大生成数量，防止爆炸
        }
        int local_total = 0;
        vector<string> local_guesses;
        #pragma omp parallel
        {
            vector<string> thread_guesses;
            int thread_total = 0;
            #pragma omp for nowait
            for (int i = 0; i < n; i++)
            {
                thread_guesses.emplace_back(a->ordered_values[i]);
                thread_total++;
            }
            #pragma omp critical
            {
                guesses.insert(guesses.end(), thread_guesses.begin(), thread_guesses.end());
                total_guesses += thread_total;
            }
        }
    }
    else
    {
        string guess;
        int seg_idx = 0;
        // 这个for循环的作用：给当前PT的所有segment赋予实际的值（最后一个segment除外）
        // segment值根据curr_indices中对应的值加以确定
        // 这个for循环你看不懂也没太大问题，并行算法不涉及这里的加速
        for (int idx : pt.curr_indices)
        {
            if (pt.content[seg_idx].type == 1)
            {
                guess += m.letters[m.FindLetter(pt.content[seg_idx])].ordered_values[idx];
            }
            if (pt.content[seg_idx].type == 2)
            {
                guess += m.digits[m.FindDigit(pt.content[seg_idx])].ordered_values[idx];
            }
            if (pt.content[seg_idx].type == 3)
            {
                guess += m.symbols[m.FindSymbol(pt.content[seg_idx])].ordered_values[idx];
            }
            seg_idx += 1;
            if (seg_idx == pt.content.size() - 1)
            {
                break;
            }
        }

        // 指向最后一个segment的指针，这个指针实际指向模型中的统计数据
        segment *a;
        if (pt.content[pt.content.size() - 1].type == 1)
        {
            a = &m.letters[m.FindLetter(pt.content[pt.content.size() - 1])];
        }
        if (pt.content[pt.content.size() - 1].type == 2)
        {
            a = &m.digits[m.FindDigit(pt.content[pt.content.size() - 1])];
        }
        if (pt.content[pt.content.size() - 1].type == 3)
        {
            a = &m.symbols[m.FindSymbol(pt.content[pt.content.size() - 1])];
        }
        
        // 优化：ordered_values已按概率降序排列，直接顺序并行遍历即可，保证高概率优先
        int n = pt.max_indices[pt.content.size() - 1];
        if (n > max_total)
        {
            n = max_total; // 限制单个PT的最大生成数量，防止爆炸
        }
        // n 就是当前PT的最大猜测数（即最后一个segment的所有value数量）
        #pragma omp parallel
        {
            vector<string> thread_guesses;
            int thread_total = 0;
            #pragma omp for schedule(static)
            for (int i = 0; i < n; i++)
            {
                // ordered_values[i] 概率更高的排在前面
                string temp = guess + a->ordered_values[i];
                thread_guesses.emplace_back(temp);
                thread_total++;
            }
            #pragma omp critical
            {
                guesses.insert(guesses.end(), thread_guesses.begin(), thread_guesses.end());
                total_guesses += thread_total;
            }
        }
    }
}

// 优化版Generate_Beam：只对最后2段做beam扩展，每段只取前top_k高概率value，beam_width较小
// void PriorityQueue::Generate_Beam(PT pt, int beam_width)
// {
//     const int MAX_TOTAL = 10000; // 单个PT最大生成数量
//     const int top_k = 64;       // 每段只取前top_k高概率value
//     CalProb(pt);

//     typedef pair<string, float> GuessProb;
//     string prefix;

//     int n_seg = pt.content.size();
//     // 1. 前n-2段直接拼接
//     int seg_idx = 0;
//     for (; seg_idx + 2 < n_seg; ++seg_idx) {
//         int idx = pt.curr_indices[seg_idx];
//         int type = pt.content[seg_idx].type;
//         if (type == 1) prefix += m.letters[m.FindLetter(pt.content[seg_idx])].ordered_values[idx];
//         else if (type == 2) prefix += m.digits[m.FindDigit(pt.content[seg_idx])].ordered_values[idx];
//         else if (type == 3) prefix += m.symbols[m.FindSymbol(pt.content[seg_idx])].ordered_values[idx];
//     }

//     // 2. 对倒数第2段做beam扩展
//     vector<GuessProb> beam;
//     beam.emplace_back(prefix, pt.prob);

//     if (n_seg >= 2) {
//         int seg2_idx = seg_idx;
//         int seg2_type = pt.content[seg2_idx].type;
//         int seg2_id = -1;
//         if (seg2_type == 1) seg2_id = m.FindLetter(pt.content[seg2_idx]);
//         else if (seg2_type == 2) seg2_id = m.FindDigit(pt.content[seg2_idx]);
//         else if (seg2_type == 3) seg2_id = m.FindSymbol(pt.content[seg2_idx]);
//         const vector<string>* values2 = nullptr;
//         const vector<int>* freqs2 = nullptr;
//         int total_freq2 = 1;
//         if (seg2_type == 1) { values2 = &m.letters[seg2_id].ordered_values; freqs2 = &m.letters[seg2_id].ordered_freqs; total_freq2 = m.letters[seg2_id].total_freq; }
//         else if (seg2_type == 2) { values2 = &m.digits[seg2_id].ordered_values; freqs2 = &m.digits[seg2_id].ordered_freqs; total_freq2 = m.digits[seg2_id].total_freq; }
//         else if (seg2_type == 3) { values2 = &m.symbols[seg2_id].ordered_values; freqs2 = &m.symbols[seg2_id].ordered_freqs; total_freq2 = m.symbols[seg2_id].total_freq; }
//         if (!values2 || values2->empty()) return;

//         int limit2 = std::min(top_k, (int)values2->size());
//         vector<GuessProb> next_beam;
//         for (const auto& bp : beam) {
//             for (int i = 0; i < limit2; ++i) {
//                 float prob = bp.second * ((*freqs2)[i] / float(total_freq2));
//                 next_beam.emplace_back(bp.first + (*values2)[i], prob);
//             }
//         }
//         // 剪枝
//         if ((int)next_beam.size() > beam_width) {
//             std::nth_element(next_beam.begin(), next_beam.begin() + beam_width - 1, next_beam.end(),
//                 [](const GuessProb& a, const GuessProb& b) { return a.second > b.second; });
//             next_beam.resize(beam_width);
//         }
//         beam.swap(next_beam);
//     }

//     // 3. 对最后一段做beam扩展
//     if (n_seg >= 1) {
//         int seg3_idx = n_seg - 1;
//         int seg3_type = pt.content[seg3_idx].type;
//         int seg3_id = -1;
//         if (seg3_type == 1) seg3_id = m.FindLetter(pt.content[seg3_idx]);
//         else if (seg3_type == 2) seg3_id = m.FindDigit(pt.content[seg3_idx]);
//         else if (seg3_type == 3) seg3_id = m.FindSymbol(pt.content[seg3_idx]);
//         const vector<string>* values3 = nullptr;
//         const vector<int>* freqs3 = nullptr;
//         int total_freq3 = 1;
//         if (seg3_type == 1) { values3 = &m.letters[seg3_id].ordered_values; freqs3 = &m.letters[seg3_id].ordered_freqs; total_freq3 = m.letters[seg3_id].total_freq; }
//         else if (seg3_type == 2) { values3 = &m.digits[seg3_id].ordered_values; freqs3 = &m.digits[seg3_id].ordered_freqs; total_freq3 = m.digits[seg3_id].total_freq; }
//         else if (seg3_type == 3) { values3 = &m.symbols[seg3_id].ordered_values; freqs3 = &m.symbols[seg3_id].ordered_freqs; total_freq3 = m.symbols[seg3_id].total_freq; }
//         if (!values3 || values3->empty()) return;

//         int limit3 = std::min(top_k, (int)values3->size());
//         vector<GuessProb> next_beam;
//         for (const auto& bp : beam) {
//             for (int i = 0; i < limit3; ++i) {
//                 float prob = bp.second * ((*freqs3)[i] / float(total_freq3));
//                 next_beam.emplace_back(bp.first + (*values3)[i], prob);
//             }
//         }
//         // 剪枝
//         if ((int)next_beam.size() > beam_width) {
//             std::nth_element(next_beam.begin(), next_beam.begin() + beam_width - 1, next_beam.end(),
//                 [](const GuessProb& a, const GuessProb& b) { return a.second > b.second; });
//             next_beam.resize(beam_width);
//         }
//         beam.swap(next_beam);
//     }

//     // 最终输出概率最高的MAX_TOTAL个猜测
//     int n = std::min((int)beam.size(), MAX_TOTAL);
//     std::sort(beam.begin(), beam.end(), [](const GuessProb& a, const GuessProb& b) { return a.second > b.second; });
//     #pragma omp parallel
//     {
//         vector<string> thread_guesses;
//         int thread_total = 0;
//         #pragma omp for nowait
//         for (int i = 0; i < n; ++i) {
//             thread_guesses.emplace_back(beam[i].first);
//             thread_total++;
//         }
//         #pragma omp critical
//         {
//             guesses.insert(guesses.end(), thread_guesses.begin(), thread_guesses.end());
//             total_guesses += thread_total;
//         }
//     }
// }



void PriorityQueue::Generate_Beam(PT pt, int beam_width)
{
    CalProb(pt);
    typedef pair<string, float> GuessProb;
    string prefix;
    int n_seg = pt.content.size();

    // ===== 动态 top_k 策略 =====
    auto GetDynamicTopK = [&](int value_size) -> int {
        if (value_size <= 10) return value_size;
        if (value_size <= 50) return std::min(beam_width, value_size);
        if (value_size <= 200) return std::min(beam_width, int(value_size * 0.5));
        if (value_size <= 500) return std::min(beam_width, int(value_size * 0.3));
        return std::min(beam_width, int(20 + 30 * log10(value_size)));
    };

    // ===== 1. 拼接前 n-2 段 =====
    int seg_idx = 0;
    for (; seg_idx + 2 < n_seg; ++seg_idx) {
        int idx = pt.curr_indices[seg_idx];
        int type = pt.content[seg_idx].type;
        if (type == 1) prefix += m.letters[m.FindLetter(pt.content[seg_idx])].ordered_values[idx];
        else if (type == 2) prefix += m.digits[m.FindDigit(pt.content[seg_idx])].ordered_values[idx];
        else if (type == 3) prefix += m.symbols[m.FindSymbol(pt.content[seg_idx])].ordered_values[idx];
    }

    // ===== 2. 倒数第 2 段 Beam 扩展 =====
    vector<GuessProb> beam = { {prefix, pt.prob} };

    if (n_seg >= 2) {
        int seg2_idx = seg_idx;
        int seg2_type = pt.content[seg2_idx].type;
        int seg2_id = (seg2_type == 1) ? m.FindLetter(pt.content[seg2_idx]) :
                       (seg2_type == 2) ? m.FindDigit(pt.content[seg2_idx]) :
                       m.FindSymbol(pt.content[seg2_idx]);

        const vector<string>& values2 = (seg2_type == 1) ? m.letters[seg2_id].ordered_values :
                                        (seg2_type == 2) ? m.digits[seg2_id].ordered_values :
                                                           m.symbols[seg2_id].ordered_values;
        const vector<int>& freqs2 = (seg2_type == 1) ? m.letters[seg2_id].ordered_freqs :
                                      (seg2_type == 2) ? m.digits[seg2_id].ordered_freqs :
                                                         m.symbols[seg2_id].ordered_freqs;
        int total_freq2 = (seg2_type == 1) ? m.letters[seg2_id].total_freq :
                          (seg2_type == 2) ? m.digits[seg2_id].total_freq :
                                             m.symbols[seg2_id].total_freq;

        if (values2.empty()) return;

        int limit2 = GetDynamicTopK(values2.size());
        vector<GuessProb> next_beam;
        for (const auto& bp : beam) {
            for (int i = 0; i < limit2; ++i) {
                float prob = bp.second * (float(freqs2[i]) / total_freq2);
                next_beam.emplace_back(bp.first + values2[i], prob);
            }
        }

        if ((int)next_beam.size() > beam_width) {
            std::nth_element(next_beam.begin(), next_beam.begin() + beam_width - 1, next_beam.end(),
                [](const GuessProb& a, const GuessProb& b) { return a.second > b.second; });
            next_beam.resize(beam_width);
        }
        beam.swap(next_beam);
    }

    // ===== 3. 最后一段 Beam 扩展 =====
    if (n_seg >= 1) {
        int seg3_idx = n_seg - 1;
        int seg3_type = pt.content[seg3_idx].type;
        int seg3_id = (seg3_type == 1) ? m.FindLetter(pt.content[seg3_idx]) :
                       (seg3_type == 2) ? m.FindDigit(pt.content[seg3_idx]) :
                       m.FindSymbol(pt.content[seg3_idx]);

        const vector<string>& values3 = (seg3_type == 1) ? m.letters[seg3_id].ordered_values :
                                        (seg3_type == 2) ? m.digits[seg3_id].ordered_values :
                                                           m.symbols[seg3_id].ordered_values;
        const vector<int>& freqs3 = (seg3_type == 1) ? m.letters[seg3_id].ordered_freqs :
                                      (seg3_type == 2) ? m.digits[seg3_id].ordered_freqs :
                                                         m.symbols[seg3_id].ordered_freqs;
        int total_freq3 = (seg3_type == 1) ? m.letters[seg3_id].total_freq :
                          (seg3_type == 2) ? m.digits[seg3_id].total_freq :
                                             m.symbols[seg3_id].total_freq;

        if (values3.empty()) return;

        int limit3 = GetDynamicTopK(values3.size());
        vector<GuessProb> next_beam;
        for (const auto& bp : beam) {
            for (int i = 0; i < limit3; ++i) {
                float prob = bp.second * (float(freqs3[i]) / total_freq3);
                next_beam.emplace_back(bp.first + values3[i], prob);
            }
        }

        if ((int)next_beam.size() > beam_width) {
            std::nth_element(next_beam.begin(), next_beam.begin() + beam_width - 1, next_beam.end(),
                [](const GuessProb& a, const GuessProb& b) { return a.second > b.second; });
            next_beam.resize(beam_width);
        }
        beam.swap(next_beam);
    }

    // ===== 4. 截断输出范围，避免概率分布歪斜 =====
    if (beam.empty()) return;

    std::sort(beam.begin(), beam.end(), [](const GuessProb& a, const GuessProb& b) {
        return a.second > b.second;
    });

    float cumulative = 0;
    int cutoff = beam.size();  // 默认保留所有
    for (int i = 0; i < (int)beam.size(); ++i) {
        cumulative += beam[i].second;
        if (cumulative >= 0.95f) {
            cutoff = i + 1;
            break;
        }
    }

    int min_limit = 20;
    int safe_limit = beam_width * 10;
    int n = std::min(std::max(cutoff, min_limit), safe_limit);
    n = std::min(n, (int)beam.size());  // 防止越界

    // ===== 5. 并行写入猜测结果 =====
    #pragma omp parallel
    {
        vector<string> thread_guesses;
        int thread_total = 0;

        #pragma omp for nowait
        for (int i = 0; i < n; ++i) {
            if (!beam[i].first.empty()) {
                thread_guesses.emplace_back(beam[i].first);
                thread_total++;
            }
        }

        #pragma omp critical
        {
            guesses.insert(guesses.end(), thread_guesses.begin(), thread_guesses.end());
            total_guesses += thread_total;
        }
    }
}
=======
#include "PCFG.h"
#include <queue>
#include <string_view>
#include <unordered_map>
#include <cmath>
#include <memory_resource>
#include <algorithm>
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
    size_t n = m.ordered_pts.size();
    vector<PT> pts(n);
    
    // 预计算 segment 索引，避免重复调用 Find* 函数
    vector<vector<int>> letter_indices(m.letters.size());
    vector<vector<int>> digit_indices(m.digits.size());
    vector<vector<int>> symbol_indices(m.symbols.size());
    
    #pragma omp parallel for
    for (size_t i = 0; i < m.letters.size(); ++i) {
        letter_indices[i].resize(m.letters[i].ordered_values.size());
        for (size_t j = 0; j < m.letters[i].ordered_values.size(); ++j) {
            letter_indices[i][j] = i; // 缓存索引以加快查找
        }
    }
    #pragma omp parallel for
    for (size_t i = 0; i < m.digits.size(); ++i) {
        digit_indices[i].resize(m.digits[i].ordered_values.size());
        for (size_t j = 0; j < m.digits[i].ordered_values.size(); ++j) {
            digit_indices[i][j] = i;
        }
    }
    #pragma omp parallel for
    for (size_t i = 0; i < m.symbols.size(); ++i) {
        symbol_indices[i].resize(m.symbols[i].ordered_values.size());
        for (size_t j = 0; j < m.symbols[i].ordered_values.size(); ++j) {
            symbol_indices[i][j] = i;
        }
    }

    // 并行初始化 PT 对象
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < n; ++i) {
        PT pt = m.ordered_pts[i];
        pt.max_indices.reserve(pt.content.size()); // 预分配 max_indices
        pt.max_indices.clear();
        
        for (const segment& seg : pt.content) {
            if (seg.type == 1) {
                int idx = m.FindLetter(seg);
                pt.max_indices.emplace_back(m.letters[idx].ordered_values.size());
            } else if (seg.type == 2) {
                int idx = m.FindDigit(seg);
                pt.max_indices.emplace_back(m.digits[idx].ordered_values.size());
            } else if (seg.type == 3) {
                int idx = m.FindSymbol(seg);
                pt.max_indices.emplace_back(m.symbols[idx].ordered_values.size());
            }
        }
        pt.preterm_prob = static_cast<float>(m.preterm_freq[m.FindPT(pt)]) / m.total_preterm;
        CalProb(pt);
        pts[i] = pt;
    }

    // 并行插入优先队列
    int n_threads = omp_get_max_threads();
    vector<vector<PT>> local_pts(n_threads);
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        local_pts[tid].reserve(n / n_threads + 1); // 每个线程预分配空间
        #pragma omp for schedule(static)
        for (size_t i = 0; i < n; ++i) {
            local_pts[tid].push_back(pts[i]);
        }
    }

    // 串行合并 local_pts 到优先队列
    for (const auto& thread_pts : local_pts) {
        for (const auto& pt : thread_pts) {
            priority.push(pt);
        }
    }
}

void PriorityQueue::PopNext()
{
    // 动态调整弹出数量，基于线程数和队列大小
    int n_threads = omp_get_max_threads();
    int pop_count = max(32, n_threads * 4); // 动态设置，至少32或线程数的4倍
    int actual_pop = min(pop_count, static_cast<int>(priority.size()));
    
    // 预分配 popped_pts 容器
    vector<PT> popped_pts;
    popped_pts.reserve(actual_pop);
    
    // 串行弹出 PT 对象（优先队列操作难以并行）
    for (int i = 0; i < actual_pop && !priority.empty(); ++i) {
        popped_pts.push_back(priority.top());
        priority.pop();
    }

    // 并行处理 Generate
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < popped_pts.size(); ++i) {
        //Generate(popped_pts[i]);
        Generate_Beam(popped_pts[i], 1024); // 使用新的 Generate_Beam 函数
    }

    // 并行生成新 PT 并收集到线程本地容器
    vector<vector<PT>> local_new_pts(n_threads);
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        local_new_pts[tid].reserve(popped_pts.size() / n_threads + 1); // 预分配本地容器
        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < popped_pts.size(); ++i) {
            vector<PT> new_pts = popped_pts[i].NewPTs();
            for (PT& pt : new_pts) {
                CalProb(pt);
                local_new_pts[tid].push_back(pt);
            }
        }
    }

    // 串行合并到优先队列
    for (const auto& thread_pts : local_new_pts) {
        for (const auto& pt : thread_pts) {
            priority.push(pt);
        }
    }
}

// 这个函数你就算看不懂，对并行算法的实现影响也不大
// 当然如果你想做一个基于多优先队列的并行算法，可能得稍微看一看了
vector<PT> PT::NewPTs()
{
    vector<PT> res;
    if (content.size() == 1)
        return res;

    int init_pivot = pivot;
    int n = curr_indices.size() - 1;
    res.reserve(n - init_pivot); // 预分配空间

    // 并行化收益有限，这里保持串行，保证顺序和低开销
    for (int i = init_pivot; i < n; ++i)
    {
        curr_indices[i] += 1;
        if (curr_indices[i] < max_indices[i])
        {
            pivot = i;
            res.emplace_back(*this); // 只在需要时拷贝
        }
        curr_indices[i] -= 1;
    }
    pivot = init_pivot;
    return res;
}


// 这个函数是PCFG并行化算法的主要载体
// 尽量看懂，然后进行并行实现
void PriorityQueue::Generate(PT pt)
{
    int max_total = 1000000; // 单个PT最大生成数量，建议根据实际需求调整
    // 计算PT的概率，这里主要是给PT的概率进行初始化
    CalProb(pt);

    // 对于只有一个segment的PT，直接遍历生成其中的所有value即可
    if (pt.content.size() == 1)
    {
        // 指向最后一个segment的指针，这个指针实际指向模型中的统计数据
        segment *a;
        // 在模型中定位到这个segment
        if (pt.content[0].type == 1)
        {
            a = &m.letters[m.FindLetter(pt.content[0])];
        }
        if (pt.content[0].type == 2)
        {
            a = &m.digits[m.FindDigit(pt.content[0])];
        }
        if (pt.content[0].type == 3)
        {
            a = &m.symbols[m.FindSymbol(pt.content[0])];
        }
        
        // Multi-thread TODO：
        // 这个for循环就是你需要进行并行化的主要部分了，特别是在多线程&GPU编程任务中
        // 可以看到，这个循环本质上就是把模型中一个segment的所有value，赋值到PT中，形成一系列新的猜测
        // 这个过程是可以高度并行化的
        int n = pt.max_indices[0];
        // n 就是当前PT的最大猜测数（即该segment的所有value数量）
        if (n > max_total)
        {
            n = max_total; // 限制单个PT的最大生成数量，防止爆炸
        }
        int local_total = 0;
        vector<string> local_guesses;
        #pragma omp parallel
        {
            vector<string> thread_guesses;
            int thread_total = 0;
            #pragma omp for nowait
            for (int i = 0; i < n; i++)
            {
                thread_guesses.emplace_back(a->ordered_values[i]);
                thread_total++;
            }
            #pragma omp critical
            {
                guesses.insert(guesses.end(), thread_guesses.begin(), thread_guesses.end());
                total_guesses += thread_total;
            }
        }
    }
    else
    {
        string guess;
        int seg_idx = 0;
        // 这个for循环的作用：给当前PT的所有segment赋予实际的值（最后一个segment除外）
        // segment值根据curr_indices中对应的值加以确定
        // 这个for循环你看不懂也没太大问题，并行算法不涉及这里的加速
        for (int idx : pt.curr_indices)
        {
            if (pt.content[seg_idx].type == 1)
            {
                guess += m.letters[m.FindLetter(pt.content[seg_idx])].ordered_values[idx];
            }
            if (pt.content[seg_idx].type == 2)
            {
                guess += m.digits[m.FindDigit(pt.content[seg_idx])].ordered_values[idx];
            }
            if (pt.content[seg_idx].type == 3)
            {
                guess += m.symbols[m.FindSymbol(pt.content[seg_idx])].ordered_values[idx];
            }
            seg_idx += 1;
            if (seg_idx == pt.content.size() - 1)
            {
                break;
            }
        }

        // 指向最后一个segment的指针，这个指针实际指向模型中的统计数据
        segment *a;
        if (pt.content[pt.content.size() - 1].type == 1)
        {
            a = &m.letters[m.FindLetter(pt.content[pt.content.size() - 1])];
        }
        if (pt.content[pt.content.size() - 1].type == 2)
        {
            a = &m.digits[m.FindDigit(pt.content[pt.content.size() - 1])];
        }
        if (pt.content[pt.content.size() - 1].type == 3)
        {
            a = &m.symbols[m.FindSymbol(pt.content[pt.content.size() - 1])];
        }
        
        // 优化：ordered_values已按概率降序排列，直接顺序并行遍历即可，保证高概率优先
        int n = pt.max_indices[pt.content.size() - 1];
        if (n > max_total)
        {
            n = max_total; // 限制单个PT的最大生成数量，防止爆炸
        }
        // n 就是当前PT的最大猜测数（即最后一个segment的所有value数量）
        #pragma omp parallel
        {
            vector<string> thread_guesses;
            int thread_total = 0;
            #pragma omp for schedule(static)
            for (int i = 0; i < n; i++)
            {
                // ordered_values[i] 概率更高的排在前面
                string temp = guess + a->ordered_values[i];
                thread_guesses.emplace_back(temp);
                thread_total++;
            }
            #pragma omp critical
            {
                guesses.insert(guesses.end(), thread_guesses.begin(), thread_guesses.end());
                total_guesses += thread_total;
            }
        }
    }
}

// 优化版Generate_Beam：只对最后2段做beam扩展，每段只取前top_k高概率value，beam_width较小
// void PriorityQueue::Generate_Beam(PT pt, int beam_width)
// {
//     const int MAX_TOTAL = 10000; // 单个PT最大生成数量
//     const int top_k = 64;       // 每段只取前top_k高概率value
//     CalProb(pt);

//     typedef pair<string, float> GuessProb;
//     string prefix;

//     int n_seg = pt.content.size();
//     // 1. 前n-2段直接拼接
//     int seg_idx = 0;
//     for (; seg_idx + 2 < n_seg; ++seg_idx) {
//         int idx = pt.curr_indices[seg_idx];
//         int type = pt.content[seg_idx].type;
//         if (type == 1) prefix += m.letters[m.FindLetter(pt.content[seg_idx])].ordered_values[idx];
//         else if (type == 2) prefix += m.digits[m.FindDigit(pt.content[seg_idx])].ordered_values[idx];
//         else if (type == 3) prefix += m.symbols[m.FindSymbol(pt.content[seg_idx])].ordered_values[idx];
//     }

//     // 2. 对倒数第2段做beam扩展
//     vector<GuessProb> beam;
//     beam.emplace_back(prefix, pt.prob);

//     if (n_seg >= 2) {
//         int seg2_idx = seg_idx;
//         int seg2_type = pt.content[seg2_idx].type;
//         int seg2_id = -1;
//         if (seg2_type == 1) seg2_id = m.FindLetter(pt.content[seg2_idx]);
//         else if (seg2_type == 2) seg2_id = m.FindDigit(pt.content[seg2_idx]);
//         else if (seg2_type == 3) seg2_id = m.FindSymbol(pt.content[seg2_idx]);
//         const vector<string>* values2 = nullptr;
//         const vector<int>* freqs2 = nullptr;
//         int total_freq2 = 1;
//         if (seg2_type == 1) { values2 = &m.letters[seg2_id].ordered_values; freqs2 = &m.letters[seg2_id].ordered_freqs; total_freq2 = m.letters[seg2_id].total_freq; }
//         else if (seg2_type == 2) { values2 = &m.digits[seg2_id].ordered_values; freqs2 = &m.digits[seg2_id].ordered_freqs; total_freq2 = m.digits[seg2_id].total_freq; }
//         else if (seg2_type == 3) { values2 = &m.symbols[seg2_id].ordered_values; freqs2 = &m.symbols[seg2_id].ordered_freqs; total_freq2 = m.symbols[seg2_id].total_freq; }
//         if (!values2 || values2->empty()) return;

//         int limit2 = std::min(top_k, (int)values2->size());
//         vector<GuessProb> next_beam;
//         for (const auto& bp : beam) {
//             for (int i = 0; i < limit2; ++i) {
//                 float prob = bp.second * ((*freqs2)[i] / float(total_freq2));
//                 next_beam.emplace_back(bp.first + (*values2)[i], prob);
//             }
//         }
//         // 剪枝
//         if ((int)next_beam.size() > beam_width) {
//             std::nth_element(next_beam.begin(), next_beam.begin() + beam_width - 1, next_beam.end(),
//                 [](const GuessProb& a, const GuessProb& b) { return a.second > b.second; });
//             next_beam.resize(beam_width);
//         }
//         beam.swap(next_beam);
//     }

//     // 3. 对最后一段做beam扩展
//     if (n_seg >= 1) {
//         int seg3_idx = n_seg - 1;
//         int seg3_type = pt.content[seg3_idx].type;
//         int seg3_id = -1;
//         if (seg3_type == 1) seg3_id = m.FindLetter(pt.content[seg3_idx]);
//         else if (seg3_type == 2) seg3_id = m.FindDigit(pt.content[seg3_idx]);
//         else if (seg3_type == 3) seg3_id = m.FindSymbol(pt.content[seg3_idx]);
//         const vector<string>* values3 = nullptr;
//         const vector<int>* freqs3 = nullptr;
//         int total_freq3 = 1;
//         if (seg3_type == 1) { values3 = &m.letters[seg3_id].ordered_values; freqs3 = &m.letters[seg3_id].ordered_freqs; total_freq3 = m.letters[seg3_id].total_freq; }
//         else if (seg3_type == 2) { values3 = &m.digits[seg3_id].ordered_values; freqs3 = &m.digits[seg3_id].ordered_freqs; total_freq3 = m.digits[seg3_id].total_freq; }
//         else if (seg3_type == 3) { values3 = &m.symbols[seg3_id].ordered_values; freqs3 = &m.symbols[seg3_id].ordered_freqs; total_freq3 = m.symbols[seg3_id].total_freq; }
//         if (!values3 || values3->empty()) return;

//         int limit3 = std::min(top_k, (int)values3->size());
//         vector<GuessProb> next_beam;
//         for (const auto& bp : beam) {
//             for (int i = 0; i < limit3; ++i) {
//                 float prob = bp.second * ((*freqs3)[i] / float(total_freq3));
//                 next_beam.emplace_back(bp.first + (*values3)[i], prob);
//             }
//         }
//         // 剪枝
//         if ((int)next_beam.size() > beam_width) {
//             std::nth_element(next_beam.begin(), next_beam.begin() + beam_width - 1, next_beam.end(),
//                 [](const GuessProb& a, const GuessProb& b) { return a.second > b.second; });
//             next_beam.resize(beam_width);
//         }
//         beam.swap(next_beam);
//     }

//     // 最终输出概率最高的MAX_TOTAL个猜测
//     int n = std::min((int)beam.size(), MAX_TOTAL);
//     std::sort(beam.begin(), beam.end(), [](const GuessProb& a, const GuessProb& b) { return a.second > b.second; });
//     #pragma omp parallel
//     {
//         vector<string> thread_guesses;
//         int thread_total = 0;
//         #pragma omp for nowait
//         for (int i = 0; i < n; ++i) {
//             thread_guesses.emplace_back(beam[i].first);
//             thread_total++;
//         }
//         #pragma omp critical
//         {
//             guesses.insert(guesses.end(), thread_guesses.begin(), thread_guesses.end());
//             total_guesses += thread_total;
//         }
//     }
// }



void PriorityQueue::Generate_Beam(PT pt, int beam_width)
{
    CalProb(pt);
    typedef pair<string, float> GuessProb;
    string prefix;
    int n_seg = pt.content.size();

    // ===== 动态 top_k 策略 =====
    auto GetDynamicTopK = [&](int value_size) -> int {
        if (value_size <= 10) return value_size;
        if (value_size <= 50) return std::min(beam_width, value_size);
        if (value_size <= 200) return std::min(beam_width, int(value_size * 0.5));
        if (value_size <= 500) return std::min(beam_width, int(value_size * 0.3));
        return std::min(beam_width, int(20 + 30 * log10(value_size)));
    };

    // ===== 1. 拼接前 n-2 段 =====
    int seg_idx = 0;
    for (; seg_idx + 2 < n_seg; ++seg_idx) {
        int idx = pt.curr_indices[seg_idx];
        int type = pt.content[seg_idx].type;
        if (type == 1) prefix += m.letters[m.FindLetter(pt.content[seg_idx])].ordered_values[idx];
        else if (type == 2) prefix += m.digits[m.FindDigit(pt.content[seg_idx])].ordered_values[idx];
        else if (type == 3) prefix += m.symbols[m.FindSymbol(pt.content[seg_idx])].ordered_values[idx];
    }

    // ===== 2. 倒数第 2 段 Beam 扩展 =====
    vector<GuessProb> beam = { {prefix, pt.prob} };

    if (n_seg >= 2) {
        int seg2_idx = seg_idx;
        int seg2_type = pt.content[seg2_idx].type;
        int seg2_id = (seg2_type == 1) ? m.FindLetter(pt.content[seg2_idx]) :
                       (seg2_type == 2) ? m.FindDigit(pt.content[seg2_idx]) :
                       m.FindSymbol(pt.content[seg2_idx]);

        const vector<string>& values2 = (seg2_type == 1) ? m.letters[seg2_id].ordered_values :
                                        (seg2_type == 2) ? m.digits[seg2_id].ordered_values :
                                                           m.symbols[seg2_id].ordered_values;
        const vector<int>& freqs2 = (seg2_type == 1) ? m.letters[seg2_id].ordered_freqs :
                                      (seg2_type == 2) ? m.digits[seg2_id].ordered_freqs :
                                                         m.symbols[seg2_id].ordered_freqs;
        int total_freq2 = (seg2_type == 1) ? m.letters[seg2_id].total_freq :
                          (seg2_type == 2) ? m.digits[seg2_id].total_freq :
                                             m.symbols[seg2_id].total_freq;

        if (values2.empty()) return;

        int limit2 = GetDynamicTopK(values2.size());
        vector<GuessProb> next_beam;
        for (const auto& bp : beam) {
            for (int i = 0; i < limit2; ++i) {
                float prob = bp.second * (float(freqs2[i]) / total_freq2);
                next_beam.emplace_back(bp.first + values2[i], prob);
            }
        }

        if ((int)next_beam.size() > beam_width) {
            std::nth_element(next_beam.begin(), next_beam.begin() + beam_width - 1, next_beam.end(),
                [](const GuessProb& a, const GuessProb& b) { return a.second > b.second; });
            next_beam.resize(beam_width);
        }
        beam.swap(next_beam);
    }

    // ===== 3. 最后一段 Beam 扩展 =====
    if (n_seg >= 1) {
        int seg3_idx = n_seg - 1;
        int seg3_type = pt.content[seg3_idx].type;
        int seg3_id = (seg3_type == 1) ? m.FindLetter(pt.content[seg3_idx]) :
                       (seg3_type == 2) ? m.FindDigit(pt.content[seg3_idx]) :
                       m.FindSymbol(pt.content[seg3_idx]);

        const vector<string>& values3 = (seg3_type == 1) ? m.letters[seg3_id].ordered_values :
                                        (seg3_type == 2) ? m.digits[seg3_id].ordered_values :
                                                           m.symbols[seg3_id].ordered_values;
        const vector<int>& freqs3 = (seg3_type == 1) ? m.letters[seg3_id].ordered_freqs :
                                      (seg3_type == 2) ? m.digits[seg3_id].ordered_freqs :
                                                         m.symbols[seg3_id].ordered_freqs;
        int total_freq3 = (seg3_type == 1) ? m.letters[seg3_id].total_freq :
                          (seg3_type == 2) ? m.digits[seg3_id].total_freq :
                                             m.symbols[seg3_id].total_freq;

        if (values3.empty()) return;

        int limit3 = GetDynamicTopK(values3.size());
        vector<GuessProb> next_beam;
        for (const auto& bp : beam) {
            for (int i = 0; i < limit3; ++i) {
                float prob = bp.second * (float(freqs3[i]) / total_freq3);
                next_beam.emplace_back(bp.first + values3[i], prob);
            }
        }

        if ((int)next_beam.size() > beam_width) {
            std::nth_element(next_beam.begin(), next_beam.begin() + beam_width - 1, next_beam.end(),
                [](const GuessProb& a, const GuessProb& b) { return a.second > b.second; });
            next_beam.resize(beam_width);
        }
        beam.swap(next_beam);
    }

    // ===== 4. 截断输出范围，避免概率分布歪斜 =====
    if (beam.empty()) return;

    std::sort(beam.begin(), beam.end(), [](const GuessProb& a, const GuessProb& b) {
        return a.second > b.second;
    });

    float cumulative = 0;
    int cutoff = beam.size();  // 默认保留所有
    for (int i = 0; i < (int)beam.size(); ++i) {
        cumulative += beam[i].second;
        if (cumulative >= 0.95f) {
            cutoff = i + 1;
            break;
        }
    }

    int min_limit = 20;
    int safe_limit = beam_width * 10;
    int n = std::min(std::max(cutoff, min_limit), safe_limit);
    n = std::min(n, (int)beam.size());  // 防止越界

    // ===== 5. 并行写入猜测结果 =====
    #pragma omp parallel
    {
        vector<string> thread_guesses;
        int thread_total = 0;

        #pragma omp for nowait
        for (int i = 0; i < n; ++i) {
            if (!beam[i].first.empty()) {
                thread_guesses.emplace_back(beam[i].first);
                thread_total++;
            }
        }

        #pragma omp critical
        {
            guesses.insert(guesses.end(), thread_guesses.begin(), thread_guesses.end());
            total_guesses += thread_total;
        }
    }
}
>>>>>>> 414c772c5f5b372157d1cd655d0951d27e7baec0
