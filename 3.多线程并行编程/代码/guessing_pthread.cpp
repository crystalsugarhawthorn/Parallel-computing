#include "PCFG.h"
#include <vector>
#include <string>
#include <iostream>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <functional> // 修复 std::function 错误
#include <algorithm> // 添加此行以修复 std::remove_if 未定义的错误
using namespace std;

// 用户可在此处直接调整最大线程数
int MAX_THREADS = 8; 

// 线程池实现
class ThreadPool {
public:
    ThreadPool(int num_threads) : stop(false) {
        for (int i = 0; i < num_threads; ++i) {
            workers.emplace_back([this]() {
                while (true) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(this->queue_mutex);
                        this->condition.wait(lock, [this]() { return this->stop || !this->tasks.empty(); });
                        if (this->stop && this->tasks.empty()) return;
                        task = std::move(this->tasks.front());
                        this->tasks.pop();
                    }
                    task();
                }
            });
        }
    }

    template<class F>
    void enqueue(F&& f) {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            tasks.emplace(std::forward<F>(f));
        }
        condition.notify_one();
    }

    ~ThreadPool() {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            stop = true;
        }
        condition.notify_all();
        for (std::thread &worker : workers) {
            if (worker.joinable()) worker.join();
        }
    }

private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    std::mutex queue_mutex;
    std::condition_variable condition;
    bool stop;
};

// 全局线程池实例
ThreadPool thread_pool(MAX_THREADS);

// 通用分块逻辑
std::vector<std::pair<int, int>> divide_work(int total_work, int num_threads) {
    std::vector<std::pair<int, int>> ranges;

    // 动态调整最小任务粒度，确保任务分配均匀
    int min_task_size = std::max(1000, total_work / (num_threads * 2)); 
    int threads = std::min(num_threads, total_work / min_task_size);
    if (threads == 0) threads = 1; // 至少分配一个线程

    int base = total_work / threads, extra = total_work % threads;
    int start = 0;
    for (int i = 0; i < threads; ++i) {
        int end = start + base + (i < extra ? 1 : 0);
        if (start < end) ranges.emplace_back(start, end);
        start = end;
    }
    return ranges;
}

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
            if (seg.type == 1) {
                pt.max_indices.emplace_back(m.letters[m.FindLetter(seg)].ordered_values.size());
            } else if (seg.type == 2) {
                pt.max_indices.emplace_back(m.digits[m.FindDigit(seg)].ordered_values.size());
            } else if (seg.type == 3) {
                pt.max_indices.emplace_back(m.symbols[m.FindSymbol(seg)].ordered_values.size());
            }
        }
        pt.preterm_prob = float(m.preterm_freq[m.FindPT(pt)]) / m.total_preterm;
        CalProb(pt);
        priority.emplace_back(pt);
    }
}

// PopNext 实现
void PriorityQueue::PopNext() {
    Generate(priority.front());
    vector<PT> new_pts = priority.front().NewPTs();
    for (PT pt : new_pts) {
        CalProb(pt);
        for (auto iter = priority.begin(); iter != priority.end(); iter++) {
            if (iter != priority.end() - 1 && iter != priority.begin()) {
                if (pt.prob <= iter->prob && pt.prob > (iter + 1)->prob) {
                    priority.emplace(iter + 1, pt);
                    break;
                }
            }
            if (iter == priority.end() - 1) {
                priority.emplace_back(pt);
                break;
            }
            if (iter == priority.begin() && iter->prob < pt.prob) {
                priority.emplace(iter, pt);
                break;
            }
        }
    }
    priority.erase(priority.begin());
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

// 优化后的 Generate 函数，使用线程池
void PriorityQueue::Generate(PT pt) {
    CalProb(pt);

    int total_iterations = (pt.content.size() == 1)
        ? pt.max_indices[0]
        : pt.max_indices[pt.content.size() - 1];
    size_t old_size = guesses.size();
    guesses.resize(old_size + total_iterations);

    if (total_iterations == 0) return;

    int threads = std::min(MAX_THREADS, total_iterations);
    auto ranges = divide_work(total_iterations, threads);

    // 确保任务分配均匀，避免过小任务
    ranges.erase(std::remove_if(ranges.begin(), ranges.end(), [](const std::pair<int, int>& range) {
        return (range.second - range.first) < 1000; // 移除任务数小于1000的范围
    }), ranges.end());

    if (pt.content.size() == 1) {
        segment *a;
        if (pt.content[0].type == 1)
            a = &m.letters[m.FindLetter(pt.content[0])];
        else if (pt.content[0].type == 2)
            a = &m.digits[m.FindDigit(pt.content[0])];
        else
            a = &m.symbols[m.FindSymbol(pt.content[0])];

        for (const auto& range : ranges) {
            thread_pool.enqueue([a, old_size, range, this] {
                // 日志输出
                // static std::mutex output_mutex;
                // {
                //     std::lock_guard<std::mutex> lock(output_mutex);
                //     std::cout << "Thread ID: " << std::this_thread::get_id()
                //               << " Range: [" << range.first << ", " << range.second << ")"
                //               << " Size: " << (range.second - range.first) << std::endl;
                // }
                for (int i = range.first; i < range.second; ++i) {
                    guesses[old_size + i] = a->ordered_values[i];
                }
            });
        }
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

        segment *last_seg;
        if (pt.content[pt.content.size() - 1].type == 1)
            last_seg = &m.letters[m.FindLetter(pt.content[pt.content.size() - 1])];
        else if (pt.content[pt.content.size() - 1].type == 2)
            last_seg = &m.digits[m.FindDigit(pt.content[pt.content.size() - 1])];
        else
            last_seg = &m.symbols[m.FindSymbol(pt.content[pt.content.size() - 1])];

        for (const auto& range : ranges) {
            thread_pool.enqueue([last_seg, prefix, old_size, range, this] {
                // 日志输出
                // static std::mutex output_mutex;
                // {
                //     std::lock_guard<std::mutex> lock(output_mutex);
                //     std::cout << "Thread ID: " << std::this_thread::get_id()
                //               << " Range: [" << range.first << ", " << range.second << ")"
                //               << " Size: " << (range.second - range.first) << std::endl;
                // }
                for (int i = range.first; i < range.second; ++i) {
                    guesses[old_size + i] = prefix + last_seg->ordered_values[i];
                }
            });
        }
    }

    total_guesses += total_iterations;
}

// PopNext 支持批量参数
void PriorityQueue::PopNextBatchParallel(int batch_pt_num) {
    if (priority.empty()) return;

    // 1. 选取一批 PT，数量为 batch_pt_num 或队列剩余数
    int actual_batch = std::min(batch_pt_num, static_cast<int>(priority.size()));
    if (actual_batch <= 0) return;
    std::vector<PT> batch_pts(priority.begin(), priority.begin() + actual_batch);

    // 2. 计算每个 PT 的猜测数和写入区间
    std::vector<int> pt_counts(actual_batch);
    for (int j = 0; j < actual_batch; ++j) {
        const PT& pt = batch_pts[j];
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

    // 4. 并行处理每个 PT
    int threads = std::min(MAX_THREADS, actual_batch);
    std::vector<std::thread> pool;
    pool.reserve(threads);

    std::vector<std::vector<PT>> thread_new_pts(threads); // 每个线程生成的新PT

    int pts_per_thread = (actual_batch + threads - 1) / threads; // 均匀分配任务
    for (int t = 0; t < threads; ++t) {
        int begin = t * pts_per_thread;
        int end = std::min(begin + pts_per_thread, actual_batch);
        if (begin >= end) continue; // 跳过空任务
        pool.emplace_back([this, &batch_pts, &offsets, old_size, begin, end, &thread_new_pts, t]() {
            for (int j = begin; j < end; ++j) {
                PT pt = batch_pts[j];
                int guess_offset = old_size + offsets[j];
                int guess_count = offsets[j + 1] - offsets[j];
                if (pt.content.size() == 1) {
                    // 单segment处理
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
                    // 多segment处理
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
                    if (pt.content[pt.content.size() - 1].type == 1)
                        a = &m.letters[m.FindLetter(pt.content[pt.content.size() - 1])];
                    else if (pt.content[pt.content.size() - 1].type == 2)
                        a = &m.digits[m.FindDigit(pt.content[pt.content.size() - 1])];
                    else
                        a = &m.symbols[m.FindSymbol(pt.content[pt.content.size() - 1])];
                    for (int i = 0; i < guess_count; ++i)
                        guesses[guess_offset + i] = prefix + a->ordered_values[i];
                }

                // 生成新PT
                std::vector<PT> new_pts = pt.NewPTs();
                for (PT& npt : new_pts) {
                    CalProb(npt);
                    thread_new_pts[t].push_back(std::move(npt));
                }
            }
        });
    }
    for (auto& th : pool) th.join();
    total_guesses += new_total;

    // 5. 主线程统一插入新PT到优先队列
    std::vector<PT> new_pts;
    for (const auto& thread_pts : thread_new_pts) {
        new_pts.insert(new_pts.end(), thread_pts.begin(), thread_pts.end());
    }
    priority.erase(priority.begin(), priority.begin() + actual_batch);
    priority.insert(priority.end(), new_pts.begin(), new_pts.end());
}
