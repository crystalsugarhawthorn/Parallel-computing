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

// mpi 编译指令如下
// mpic++ correctness_guess_mpi.cpp train.cpp guessing_mpi.cpp md5.cpp -o main -O2
// mpic++ correctness_guess_mpi.cpp train.cpp guessing_mpi.cpp md5.cpp -o main -O0
// mpiexec -np 8 ./main
// qsub qsub_mpi.sh

// 最终修正版：通过发送每个向量的显式大小，使广播协议更健壮，能处理数据不一致的情况。
void BroadcastSegmentsOptimized(vector<segment> &segments, int root_rank) {
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    vector<int> int_buffer;
    vector<char> char_buffer;

    // --- 1. 序列化 (仅在 root_rank 执行) ---
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

    // --- 2. 通信 (此部分无变化) ---
    long long int_buffer_size = (world_rank == root_rank) ? int_buffer.size() : 0;
    MPI_Bcast(&int_buffer_size, 1, MPI_LONG_LONG, root_rank, MPI_COMM_WORLD);
    if (world_rank != root_rank) int_buffer.resize(int_buffer_size);
    MPI_Bcast(int_buffer.data(), int_buffer_size, MPI_INT, root_rank, MPI_COMM_WORLD);

    long long char_buffer_size = (world_rank == root_rank) ? char_buffer.size() : 0;
    MPI_Bcast(&char_buffer_size, 1, MPI_LONG_LONG, root_rank, MPI_COMM_WORLD);
    if (world_rank != root_rank) char_buffer.resize(char_buffer_size);
    MPI_Bcast(char_buffer.data(), char_buffer_size, MPI_CHAR, root_rank, MPI_COMM_WORLD);
    
    // --- 3. 反序列化 (在其他 rank 执行) ---
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
                        // 如果发生这种情况，说明数据传输有严重问题，打印错误并中止
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

// 其余函数保持不变
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

void BroadcastTestSet(unordered_set<string>& test_set, int root_rank) {
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    long long buffer_size = 0;
    vector<char> buffer;
    if (world_rank == root_rank) {
        for(const auto& pw : test_set) {
            buffer_size += pw.length() + 1;
        }
        buffer.resize(buffer_size);
        char* ptr = buffer.data();
        for(const auto& pw : test_set) {
            strcpy(ptr, pw.c_str());
            ptr += pw.length() + 1;
        }
    }
    MPI_Bcast(&buffer_size, 1, MPI_LONG_LONG, root_rank, MPI_COMM_WORLD);
    if(world_rank != root_rank) buffer.resize(buffer_size);
    MPI_Bcast(buffer.data(), buffer_size, MPI_CHAR, root_rank, MPI_COMM_WORLD);
    if(world_rank != root_rank) {
        test_set.clear();
        const char* ptr = buffer.data();
        const char* end_ptr = buffer.data() + buffer_size;
        while(ptr < end_ptr) {
            string pw(ptr);
            test_set.insert(pw);
            ptr += pw.length() + 1;
        }
    }
}

// main 函数完全不变
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    PriorityQueue q;
    double time_train = 0, time_hash = 0, time_cracked_check = 0;

    if (world_rank == 0) {
        double start = MPI_Wtime();
        q.m.train("/guessdata/Rockyou-singleLined-full.txt");
        q.m.order();
        time_train = MPI_Wtime() - start;
        cout << "Training completed in " << time_train << " s." << endl;
    }

    BroadcastModel(q, 0);

    unordered_set<string> test_set;
    if (world_rank == 0) {
        ifstream in("/guessdata/Rockyou-singleLined-full.txt");
        string pw;
        int count = 0;
        while (in >> pw && count++ < 1000000) test_set.insert(pw);
        cout << "Loaded " << test_set.size() << " test passwords." << endl;
    }
    BroadcastTestSet(test_set, 0);

    if (world_rank == 0) q.init();

    const int batch_pt_num = 10;
    int64_t global_total_generated = 0;
    int total_cracked = 0;
    double t_global_start = MPI_Wtime();

    while (true) {
        int keep_going = 1;
        vector<PT> batch_pts;

        if (world_rank == 0) {
            for (int i = 0; i < batch_pt_num && !q.priority.empty(); ++i) {
                batch_pts.push_back(q.priority.top());
                q.PopNext_MPI();
            }
            keep_going = !batch_pts.empty();
        }

        MPI_Bcast(&keep_going, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (!keep_going) break;

        int pt_count = batch_pts.size();
        MPI_Bcast(&pt_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (world_rank != 0) batch_pts.resize(pt_count);

        for (int i = 0; i < pt_count; ++i)
            BroadcastPT(batch_pts[i], 0);

        vector<string> local_guesses;
        for (const auto& pt : batch_pts) {
            vector<string> guesses_this_pt;
            q.Generate(pt, guesses_this_pt);
            local_guesses.insert(local_guesses.end(), guesses_this_pt.begin(), guesses_this_pt.end());
        }


        double t_hash0 = MPI_Wtime();
        for (const auto& guess : local_guesses) {
            bit32 state[4];
            MD5Hash(guess, state);
        }
        double hash_time = MPI_Wtime() - t_hash0;

        int cracked_local = 0;
        double t_check0 = MPI_Wtime();
        for (const auto& guess : local_guesses)
            if (test_set.count(guess)) cracked_local++;
        double check_time = MPI_Wtime() - t_check0;

        int64_t local_gen = local_guesses.size(), global_gen = 0;
        int cracked_total = 0;
        double max_hash = 0, max_check = 0;

        MPI_Reduce(&local_gen, &global_gen, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&cracked_local, &cracked_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&hash_time, &max_hash, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&check_time, &max_check, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        if (world_rank == 0) {
            global_total_generated += global_gen;
            total_cracked += cracked_total;
            time_hash += max_hash;
            time_cracked_check += max_check;

            if (global_total_generated % 10000 < batch_pt_num * 100) {
                cout << "Guessed: " << global_total_generated
                     << ", Cracked: " << total_cracked << endl;
            }
        }

        MPI_Bcast(&global_total_generated, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
        if (global_total_generated >= 10000000) break;
    }

    double t_global_end = MPI_Wtime();

    if (world_rank == 0) {
        double t_total = t_global_end - t_global_start;
        cout << "\n--- Report ---\n";
        cout << "Total generated: " << global_total_generated << endl;
        cout << "Total cracked:  " << total_cracked << endl;
        cout << "Train time:     " << time_train << " s\n";
        cout << "Guess time:     " << t_total << " s\n";
        cout << "Hashing time:   " << time_hash << " s\n";
        cout << "Check time:     " << time_cracked_check << " s\n";
    }

    MPI_Finalize();
    return 0;
}
