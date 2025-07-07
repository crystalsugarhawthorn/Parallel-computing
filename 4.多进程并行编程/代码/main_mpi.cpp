#include "PCFG.h"
#include "md5.h"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include <unordered_set>
#include <vector>

using namespace std;
using namespace chrono;

// MPI 编译指令
// mpic++ main_mpi.cpp train.cpp guessing_mpi.cpp md5.cpp -o main -O2 -fopenmp
// mpiexec -np 8 ./main
// qsub qsub_mpi.sh

// 使用第二个代码的优化广播方法
void BroadcastSegmentsOptimized(vector<segment> &segments, int root_rank) {
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    vector<int> int_buffer;
    vector<char> char_buffer;

    // 序列化 (仅在 root_rank 执行)
    if (world_rank == root_rank) {
        int_buffer.push_back(segments.size());

        for (const auto& seg : segments) {
            int_buffer.push_back(seg.type);
            int_buffer.push_back(seg.length);
            int_buffer.push_back(seg.total_freq);
            
            // 健壮性修正: 显式发送两个向量的大小
            int_buffer.push_back(seg.ordered_values.size());
            int_buffer.push_back(seg.ordered_freqs.size());

            for (const auto& val : seg.ordered_values) {
                int_buffer.push_back(val.length());
                char_buffer.insert(char_buffer.end(), val.begin(), val.end());
            }
            int_buffer.insert(int_buffer.end(), seg.ordered_freqs.begin(), seg.ordered_freqs.end());
        }
    }

    // 通信
    long long int_buffer_size = (world_rank == root_rank) ? int_buffer.size() : 0;
    MPI_Bcast(&int_buffer_size, 1, MPI_LONG_LONG, root_rank, MPI_COMM_WORLD);
    if (world_rank != root_rank) int_buffer.resize(int_buffer_size);
    MPI_Bcast(int_buffer.data(), int_buffer_size, MPI_INT, root_rank, MPI_COMM_WORLD);

    long long char_buffer_size = (world_rank == root_rank) ? char_buffer.size() : 0;
    MPI_Bcast(&char_buffer_size, 1, MPI_LONG_LONG, root_rank, MPI_COMM_WORLD);
    if (world_rank != root_rank) char_buffer.resize(char_buffer_size);
    MPI_Bcast(char_buffer.data(), char_buffer_size, MPI_CHAR, root_rank, MPI_COMM_WORLD);
    
    // 反序列化 (在其他 rank 执行)
    if (world_rank != root_rank) {
        segments.clear();
        if (int_buffer.empty()) return;

        size_t int_idx = 0;
        size_t char_idx = 0;

        int num_segments = int_buffer[int_idx++];
        segments.resize(num_segments);

        for (int i = 0; i < num_segments; ++i) {
            segments[i].type = int_buffer[int_idx++];
            segments[i].length = int_buffer[int_idx++];
            segments[i].total_freq = int_buffer[int_idx++];
            
            // 健壮性修正: 读取两个向量各自的大小
            int val_count = int_buffer[int_idx++];
            int freq_count = int_buffer[int_idx++];

            segments[i].ordered_values.resize(val_count);
            segments[i].ordered_freqs.resize(freq_count);

            for (int j = 0; j < val_count; ++j) {
                int len = int_buffer[int_idx++];
                if (len > 0) {
                    // 增加边界检查，防止 char_idx 越界
                    if (char_idx + len > char_buffer.size()) {
                        cerr << "Rank " << world_rank << ": Fatal error during deserialization. Char buffer overflow." << endl;
                        MPI_Abort(MPI_COMM_WORLD, -1);
                    }
                    segments[i].ordered_values[j].assign(&char_buffer[char_idx], len);
                    char_idx += len;
                }
            }
            // 使用独立的 freq_count 进行循环
            for (int j = 0; j < freq_count; ++j) {
                segments[i].ordered_freqs[j] = int_buffer[int_idx++];
            }
        }
    }
}

void BroadcastModel(PriorityQueue &q, int root_rank) {
    BroadcastSegmentsOptimized(q.m.letters, root_rank);
    BroadcastSegmentsOptimized(q.m.digits, root_rank);
    BroadcastSegmentsOptimized(q.m.symbols, root_rank);
}

void BroadcastPT(PT &pt, int root_rank) {
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    int num_seg = (world_rank == root_rank) ? pt.content.size() : 0;
    MPI_Bcast(&num_seg, 1, MPI_INT, root_rank, MPI_COMM_WORLD);
    if (world_rank != root_rank) pt.content.resize(num_seg);
    
    for (int i = 0; i < num_seg; ++i) {
        int type = 0, length = 0;
        if (world_rank == root_rank) {
            type = pt.content[i].type;
            length = pt.content[i].length;
        }
        MPI_Bcast(&type, 1, MPI_INT, root_rank, MPI_COMM_WORLD);
        MPI_Bcast(&length, 1, MPI_INT, root_rank, MPI_COMM_WORLD);
        if (world_rank != root_rank) {
            pt.content[i].type = type;
            pt.content[i].length = length;
        }
    }
    
    auto broadcast_vec = [&](vector<int> &vec) {
        int size = (world_rank == root_rank) ? vec.size() : 0;
        MPI_Bcast(&size, 1, MPI_INT, root_rank, MPI_COMM_WORLD);
        if (world_rank != root_rank) vec.resize(size);
        if (size > 0) MPI_Bcast(vec.data(), size, MPI_INT, root_rank, MPI_COMM_WORLD);
    };
    
    broadcast_vec(pt.curr_indices);
    broadcast_vec(pt.max_indices);
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    double time_guess = 0, time_train = 0, time_hash = 0;
    PriorityQueue q;

    // 训练阶段 - 只在 rank 0 执行
    if (world_rank == 0) {
        double start_train = MPI_Wtime();
        q.m.train("/guessdata/Rockyou-singleLined-full.txt");
        q.m.order();
        double end_train = MPI_Wtime();
        time_train = end_train - start_train;
        cout << "Training completed in " << time_train << " seconds." << endl;
    }

    // 广播模型到所有进程
    MPI_Barrier(MPI_COMM_WORLD);
    BroadcastModel(q, 0);
    MPI_Barrier(MPI_COMM_WORLD);

    // 保留第一个代码的测试集加载方式（仅用于计数）
    unordered_set<string> test_set;
    if (world_rank == 0) {
        ifstream test_data("/guessdata/Rockyou-singleLined-full.txt");
        string pw;
        int test_count = 0;
        while (test_data >> pw) {
            test_set.insert(pw);
            if (++test_count >= 1000000) break;
        }
        cout << "Loaded " << test_set.size() << " test passwords for reference." << endl;
    }

    // 初始化优先队列
    if (world_rank == 0) {
        q.init();
        cout << "Priority queue initialized." << endl;
    }

    // 使用第二个代码的批处理方式
    const int batch_pt_num = 10;  // 批处理大小
    int64_t global_total_generated = 0;
    int64_t print_threshold = 0;
    int64_t max_total_guesses = 10000000;

    double t_global_start = MPI_Wtime();

    while (true) {
        int keep_going = 1;
        vector<PT> batch_pts;

        // rank 0 准备一批 PT
        if (world_rank == 0) {
            for (int i = 0; i < batch_pt_num && !q.priority.empty(); ++i) {
                batch_pts.push_back(q.priority.top());
                q.PopNext_MPI();
            }
            keep_going = !batch_pts.empty();
        }

        // 广播是否继续
        MPI_Bcast(&keep_going, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (!keep_going) break;

        // 广播批次中的 PT 数量
        int pt_count = batch_pts.size();
        MPI_Bcast(&pt_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (world_rank != 0) batch_pts.resize(pt_count);

        // 广播每个 PT
        for (int i = 0; i < pt_count; ++i) {
            BroadcastPT(batch_pts[i], 0);
        }

        // 所有进程生成密码
        vector<string> local_guesses;
        for (const auto& pt : batch_pts) {
            vector<string> guesses_this_pt;
            q.Generate(pt, guesses_this_pt);
            local_guesses.insert(local_guesses.end(), guesses_this_pt.begin(), guesses_this_pt.end());
        }

        // 保留第一个代码的MD5哈希计算方式（批量处理）
        double local_hash_time = 0;
        double t_hash0 = MPI_Wtime();
        
        // 使用4个一组的批量哈希计算
        for (size_t i = 0; i + 3 < local_guesses.size(); i += 4) {
            bit32 states[4][4];
            MD5HashFour(local_guesses[i], local_guesses[i + 1], 
                       local_guesses[i + 2], local_guesses[i + 3], states);
        }
        // 处理剩余的密码
        for (size_t i = (local_guesses.size() / 4) * 4; i < local_guesses.size(); ++i) {
            bit32 state[4];
            MD5Hash(local_guesses[i], state);
        }
        
        double t_hash1 = MPI_Wtime();
        local_hash_time = t_hash1 - t_hash0;

        // 统计结果
        int64_t local_generated = local_guesses.size();
        int64_t global_generated = 0;
        double global_hash_time = 0;

        MPI_Reduce(&local_generated, &global_generated, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&local_hash_time, &global_hash_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        // 检查是否需要终止
        int terminate = 0;
        if (world_rank == 0) {
            global_total_generated += global_generated;
            time_hash += global_hash_time;

            // 保留第一个代码的打印方式
            if (global_total_generated >= print_threshold + 100000) {
                cout << "Guesses generated: " << global_total_generated << endl;
                print_threshold = global_total_generated;
            }

            if (global_total_generated >= max_total_guesses) {
                terminate = 1;
            }
        }

        MPI_Bcast(&terminate, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (terminate) break;

        MPI_Barrier(MPI_COMM_WORLD);
    }

    // 输出最终结果（保留第一个代码的输出格式）
    if (world_rank == 0) {
        double t_global_end = MPI_Wtime();
        time_guess = t_global_end - t_global_start;

        cout << "\n--- Final Report ---" << endl;
        cout << "Total generated: " << global_total_generated << endl;
        cout << "Guess time: " << time_guess - time_hash << " seconds" << endl;
        cout << "Hash time: " << time_hash << " seconds" << endl;
        cout << "Train time: " << time_train << " seconds" << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}