#include "PCFG.h"
#include <vector>
#include <string>
#include <iostream>
#include <queue>
#include <algorithm> 
#include <omp.h>     
using namespace std;

// CalProb 实现
void PriorityQueue::CalProb(PT &pt) {
    pt.prob = pt.preterm_prob;
    int index = 0;
    for (int idx : pt.curr_indices) {
        if (pt.content[index].type == 1) {
            pt.prob *= m.letters[m.FindLetter(pt.content[index])].ordered_freqs[idx];
            pt.prob /= m.letters[m.FindLetter(pt.content[index])].total_freq;
        } else if (pt.content[index].type == 2) {
            pt.prob *= m.digits[m.FindDigit(pt.content[index])].ordered_freqs[idx];
            pt.prob /= m.digits[m.FindDigit(pt.content[index])].total_freq;
        } else if (pt.content[index].type == 3) {
            pt.prob *= m.symbols[m.FindSymbol(pt.content[index])].ordered_freqs[idx];
            pt.prob /= m.symbols[m.FindSymbol(pt.content[index])].total_freq;
        }
        index += 1;
    }
}

// init 实现
void PriorityQueue::init() {
    for (PT pt : m.ordered_pts) {
        for (segment seg : pt.content) {
            if (seg.type == 1)
                pt.max_indices.emplace_back(m.letters[m.FindLetter(seg)].ordered_values.size());
            else if (seg.type == 2)
                pt.max_indices.emplace_back(m.digits[m.FindDigit(seg)].ordered_values.size());
            else if (seg.type == 3)
                pt.max_indices.emplace_back(m.symbols[m.FindSymbol(seg)].ordered_values.size());
        }
        pt.preterm_prob = float(m.preterm_freq[m.FindPT(pt)]) / m.total_preterm;
        CalProb(pt);
        priority.push(pt); // 修复 emplace_back 错误
    }
}

// 修改后的 PopNext 实现，基于优先队列（beam-search）机制
void PriorityQueue::PopNext() {
    if (priority.empty()) return;
    PT topPt = priority.top();
    priority.pop();
    Generate(topPt);
    vector<PT> new_pts = topPt.NewPTs();
    for (PT &pt : new_pts) {
        CalProb(pt);
        priority.push(pt);
    }
}

// PT::NewPTs 实现
vector<PT> PT::NewPTs() {
    vector<PT> res;
    if (content.size() == 1) {
        return res;
    } else {
        int init_pivot = pivot;
        for (int i = pivot; i < curr_indices.size() - 1; i += 1) {
            curr_indices[i] += 1;
            if (curr_indices[i] < max_indices[i]) {
                pivot = i;
                res.emplace_back(*this);
            }
            curr_indices[i] -= 1;
        }
        pivot = init_pivot;
        return res;
    }
}

void PriorityQueue::Generate(PT pt) {
    CalProb(pt);
    const int threshold = 50;
    const int n_threads = omp_get_max_threads();

    if (pt.content.size() == 1) {
        segment *a;
        if (pt.content[0].type == 1)
            a = &m.letters[m.FindLetter(pt.content[0])];
        else if (pt.content[0].type == 2)
            a = &m.digits[m.FindDigit(pt.content[0])];
        else
            a = &m.symbols[m.FindSymbol(pt.content[0])];

        const int n = pt.max_indices[0];
        const size_t old_size = guesses.size();
        guesses.resize(old_size + n);

        if (n >= threshold) {
            #pragma omp parallel for schedule(static)
            for (int i = 0; i < n; ++i) {
                guesses[old_size + i] = a->ordered_values[i];
            }
        } else {
            for (int i = 0; i < n; ++i) {
                guesses[old_size + i] = a->ordered_values[i];
            }
        }
        total_guesses += n;
    } else {
        string guess;
        int seg_idx = 0;
        for (int idx : pt.curr_indices) {
            if (pt.content[seg_idx].type == 1)
                guess += m.letters[m.FindLetter(pt.content[seg_idx])].ordered_values[idx];
            else if (pt.content[seg_idx].type == 2)
                guess += m.digits[m.FindDigit(pt.content[seg_idx])].ordered_values[idx];
            else
                guess += m.symbols[m.FindSymbol(pt.content[seg_idx])].ordered_values[idx];
            seg_idx += 1;
            if (seg_idx == pt.content.size() - 1)
                break;
        }

        segment *a;
        if (pt.content[pt.content.size() - 1].type == 1)
            a = &m.letters[m.FindLetter(pt.content[pt.content.size() - 1])];
        else if (pt.content[pt.content.size() - 1].type == 2)
            a = &m.digits[m.FindDigit(pt.content[pt.content.size() - 1])];
        else
            a = &m.symbols[m.FindSymbol(pt.content[pt.content.size() - 1])];

        const int n = pt.max_indices[pt.content.size() - 1];
        const size_t old_size = guesses.size();
        guesses.resize(old_size + n);

        if (n >= threshold) {
            #pragma omp parallel for schedule(static)
            for (int i = 0; i < n; ++i) {
                guesses[old_size + i] = guess + a->ordered_values[i];
            }
        } else {
            for (int i = 0; i < n; ++i) {
                guesses[old_size + i] = guess + a->ordered_values[i];
            }
        }
        total_guesses += n;
    }
}

// PopNext 支持批量参数
void PriorityQueue::PopNextBatchParallel(int batch_pt_num) {
    if (priority.empty()) return;

    // 1. 批量 pop 出 pt
    int actual_batch = std::min(batch_pt_num, static_cast<int>(priority.size()));
    if (actual_batch <= 0) return;
    vector<PT> batch_pts;
    for (int i = 0; i < actual_batch; ++i) {
        batch_pts.push_back(priority.top());
        priority.pop();
    }

    // 2. 计算每个 pt 的猜测数和写入区间
    std::vector<int> pt_counts(actual_batch);
    for (int j = 0; j < actual_batch; ++j) {
        const PT &pt = batch_pts[j];
        pt_counts[j] = (pt.content.size() == 1) ? pt.max_indices[0] : pt.max_indices.back();
    }
    std::vector<int> offsets(actual_batch + 1, 0);
    for (int j = 0; j < actual_batch; ++j) {
        offsets[j + 1] = offsets[j] + pt_counts[j];
    }
    int new_total = offsets.back();

    // 3. 预分配 guesses 空间
    int old_size = guesses.size();
    guesses.resize(old_size + new_total);

    // 4. 并行处理每个 pt
    int threads = omp_get_max_threads(); // 使用 OpenMP 获取最大线程数
    std::vector<std::vector<PT>> thread_new_pts(threads); // 每个线程生成的新 pt

    #pragma omp parallel for num_threads(threads)
    for (int t = 0; t < threads; ++t) {
        int begin = t * (actual_batch / threads);
        int end = (t == threads - 1) ? actual_batch : begin + (actual_batch / threads);
        for (int j = begin; j < end; ++j) {
            PT pt = batch_pts[j];
            int guess_offset = old_size + offsets[j];
            int guess_count = offsets[j + 1] - offsets[j];
            if (pt.content.size() == 1) {
                segment *a;
                if (pt.content[0].type == 1)
                    a = &m.letters[m.FindLetter(pt.content[0])];
                else if (pt.content[0].type == 2)
                    a = &m.digits[m.FindDigit(pt.content[0])];
                else
                    a = &m.symbols[m.FindSymbol(pt.content[0])];
                for (int i = 0; i < guess_count; ++i)
                    guesses[guess_offset + i] = a->ordered_values[i];
            } else {
                string prefix;
                int seg_idx = 0;
                for (int idx : pt.curr_indices) {
                    if (pt.content[seg_idx].type == 1)
                        prefix += m.letters[m.FindLetter(pt.content[seg_idx])].ordered_values[idx];
                    else if (pt.content[seg_idx].type == 2)
                        prefix += m.digits[m.FindDigit(pt.content[seg_idx])].ordered_values[idx];
                    else
                        prefix += m.symbols[m.FindSymbol(pt.content[seg_idx])].ordered_values[idx];
                    seg_idx++;
                    if (seg_idx == pt.content.size() - 1)
                        break;
                }
                segment *a;
                if (pt.content.back().type == 1)
                    a = &m.letters[m.FindLetter(pt.content.back())];
                else if (pt.content.back().type == 2)
                    a = &m.digits[m.FindDigit(pt.content.back())];
                else
                    a = &m.symbols[m.FindSymbol(pt.content.back())];
                for (int i = 0; i < guess_count; ++i)
                    guesses[guess_offset + i] = prefix + a->ordered_values[i];
            }

            // 生成新 pt
            std::vector<PT> new_pts = pt.NewPTs();
            for (PT &npt : new_pts) {
                CalProb(npt);
                #pragma omp critical
                thread_new_pts[t].push_back(std::move(npt));
            }
        }
    }

    // 5. 将所有线程产生的新 pt 统一 push 入优先队列
    vector<PT> new_pts;
    for (const auto &pts : thread_new_pts) {
        new_pts.insert(new_pts.end(), pts.begin(), pts.end());
    }
    for (const auto &pt : new_pts) {
        priority.push(pt);
    }
}
