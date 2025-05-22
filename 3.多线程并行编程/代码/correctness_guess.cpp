#include "PCFG.h"
#include <chrono>
#include <fstream>
#include "md5.h"
#include <iomanip>
#include <unordered_set>
using namespace std;
using namespace chrono;

// 编译指令如下
// g++ correctness_guess.cpp train.cpp guessing.cpp md5.cpp -o main
// g++ correctness_guess.cpp train.cpp guessing.cpp md5.cpp -o main -O1
// g++ correctness_guess.cpp train.cpp guessing.cpp md5.cpp -o main -O2
// g++ correctness_guess.cpp train.cpp guessing_openmp.cpp md5.cpp -o main -fopenmp -O2
// g++ correctness_guess.cpp train.cpp guessing_pthread.cpp md5.cpp -o main -O2
int main()
{
    double time_hash = 0;  // 用于MD5哈希的时间
    double time_guess = 0; // 哈希和猜测的总时长
    double time_train = 0; // 模型训练的总时长
    double time_cracked_check = 0; // 猜测检查的总时长
    PriorityQueue q;
    auto start_train = system_clock::now();
    q.m.train("/guessdata/Rockyou-singleLined-full.txt");
    q.m.order();
    auto end_train = system_clock::now();
    auto duration_train = duration_cast<microseconds>(end_train - start_train);
    time_train = double(duration_train.count()) * microseconds::period::num / microseconds::period::den;


    
    // 加载一些测试数据
    unordered_set<std::string> test_set;
    ifstream test_data("/guessdata/Rockyou-singleLined-full.txt");
    int test_count=0;
    string pw;
    while(test_data>>pw)
    {   
        test_count+=1;
        test_set.insert(pw);
        if (test_count>=1000000)
        {
            break;
        }
    }
    int cracked=0;

    q.init();
    cout << "here" << endl;
    int curr_num = 0;
    auto start = system_clock::now();
    // 由于需要定期清空内存，我们在这里记录已生成的猜测总数
    int history = 0;
    // std::ofstream a("./files/results.txt");
    while (!q.priority.empty())
    {
        q.PopNext();
        //q.PopNextBatchParallel(100);
        q.total_guesses = q.guesses.size();
        if (q.total_guesses - curr_num >= 100000)
        {
            cout << "Guesses generated: " <<history + q.total_guesses << endl;
            curr_num = q.total_guesses;

            // 在此处更改实验生成的猜测上限
            int generate_n=10000000;
            if (history + q.total_guesses > 10000000)
            {
                auto end = system_clock::now();
                auto duration = duration_cast<microseconds>(end - start);
                time_guess = double(duration.count()) * microseconds::period::num / microseconds::period::den;
                cout << "Guess time:" << time_guess - time_hash - time_cracked_check << " seconds" << endl;
                cout << "Hash time:" << time_hash << " seconds" << endl;
                cout << "Cracked check time:" << time_cracked_check << " seconds" << endl;
                cout << "Train time:" << time_train << " seconds" << endl;
                cout << "Cracked: " << cracked << endl;
                break;
            }
        }
        // 为了避免内存超限，我们在q.guesses中口令达到一定数目时，将其中的所有口令取出并且进行哈希
        // 然后，q.guesses将会被清空。为了有效记录已经生成的口令总数，维护一个history变量来进行记录
        if (curr_num > 1000000)
        {
            size_t i;

            // Step 1: 判断猜中（不计入 hash_time）
            // Step 1: 判断猜中（记录 cracked 时间）
            auto start_check = system_clock::now();
            for (i = 0; i + 3 < q.guesses.size(); i += 4)
            {
                cracked += test_set.count(q.guesses[i]) + test_set.count(q.guesses[i + 1]) +
                        test_set.count(q.guesses[i + 2]) + test_set.count(q.guesses[i + 3]);
            }
            for (; i < q.guesses.size(); i++)
            {
                cracked += test_set.count(q.guesses[i]);
            }
            auto end_check = system_clock::now();
            auto duration_check = duration_cast<microseconds>(end_check - start_check);
            time_cracked_check += double(duration_check.count()) * microseconds::period::num / microseconds::period::den;

            // Step 2: 单独计时 MD5（只算哈希时间）
            auto start_hash = system_clock::now();

            for (i = 0; i + 3 < q.guesses.size(); i += 4)
            {
                bit32 states[4][4];
                MD5HashFour(q.guesses[i], q.guesses[i + 1], q.guesses[i + 2], q.guesses[i + 3], states);
            }
            for (; i < q.guesses.size(); i++)
            {
                bit32 state[4];
                MD5Hash(q.guesses[i], state);
            }

            auto end_hash = system_clock::now();
            auto duration = duration_cast<microseconds>(end_hash - start_hash);
            time_hash += double(duration.count()) * microseconds::period::num / microseconds::period::den;

            // 记录状态
            history += curr_num;
            curr_num = 0;
            q.guesses.clear();
        }
    }
}
