<<<<<<< HEAD
#include "PCFG.h"
#include <fstream>
#include <cctype>
#include <algorithm>
#include <omp.h>
#include <vector>
#include <unordered_map>

// 这个文件里面的各函数你都不需要完全理解，甚至根本不需要看
// 从学术价值上讲，加速模型的训练过程是一个没什么价值的问题，因为我们一般假定统计学模型的训练成本较低
// 但是，假如你是一个投稿时顶着ddl做实验的倒霉研究生/实习生，提高训练速度就可以大幅节省你的时间了
// 所以如果你愿意，也可以尝试用多线程加速训练过程

/**
 * 怎么加速PCFG训练过程？据助教所知，没有公开文献提出过有效的加速方法（因为这么做基本无学术价值）
 * 
 * 但是统计学模型好就好在其数据是可加的。例如，假如我把数据集拆分成4个部分，并行训练4个不同的模型。
 * 然后我可以直接将四个模型的统计数据进行简单加和，就得到了和串行训练相同的模型了。
 * 
 * 说起来容易，做起来不一定容易，你可能会碰到一系列具体的工程问题。如果你决定加速训练过程，祝你好运！
 * 
 */

// 训练的wrapper，实际上就是读取训练集
void model::train(string path)
{
    string pw;
    ifstream train_set(path);
    int lines = 0;
    cout<<"Training..."<<endl;
    cout<<"Training phase 1: reading and parsing passwords..."<<endl;

    // 1. 先把所有密码读入内存
    vector<string> passwords;
    passwords.reserve(3000000);
    while (train_set >> pw)
    {
        lines += 1;
        if (lines % 100000 == 0)
        {
            cout <<"Lines processed: "<< lines << endl;
            if (lines > 3000000)
            {
                break;
            }
        }
        passwords.push_back(pw);
    }

    int n_threads = omp_get_max_threads();
    vector<model> local_models(n_threads);

    // 2. 流水线并行处理：分块读取、分块解析、分块合并
    const size_t block_size = 100000; // 可根据内存和核数调整
    size_t total = passwords.size();
    size_t num_blocks = (total + block_size - 1) / block_size;

    // 每个线程有自己的本地桶，减少同步
    vector<vector<PT>> pt_blocks(num_blocks);
    vector<vector<int>> pt_freq_blocks(num_blocks);
    vector<vector<segment>> letter_blocks(num_blocks), digit_blocks(num_blocks), symbol_blocks(num_blocks);
    vector<vector<int>> letter_freq_blocks(num_blocks), digit_freq_blocks(num_blocks), symbol_freq_blocks(num_blocks);
    vector<int> total_preterm_blocks(num_blocks, 0);

    #pragma omp parallel for schedule(dynamic)
    for (size_t b = 0; b < num_blocks; ++b) {
        size_t start = b * block_size;
        size_t end = min(start + block_size, total);
        model local_model;
        local_model.preterminals.reserve(block_size / 2);
        local_model.letters.reserve(128);
        local_model.digits.reserve(128);
        local_model.symbols.reserve(128);
        for (size_t i = start; i < end; ++i) {
            local_model.parse(passwords[i]);
        }
        // 收集本块数据
        pt_blocks[b] = local_model.preterminals;
        pt_freq_blocks[b].resize(local_model.preterminals.size());
        for (size_t i = 0; i < local_model.preterminals.size(); ++i)
            pt_freq_blocks[b][i] = local_model.preterm_freq[i];
        letter_blocks[b] = local_model.letters;
        letter_freq_blocks[b].resize(local_model.letters.size());
        for (size_t i = 0; i < local_model.letters.size(); ++i)
            letter_freq_blocks[b][i] = local_model.letters_freq[i];
        digit_blocks[b] = local_model.digits;
        digit_freq_blocks[b].resize(local_model.digits.size());
        for (size_t i = 0; i < local_model.digits.size(); ++i)
            digit_freq_blocks[b][i] = local_model.digits_freq[i];
        symbol_blocks[b] = local_model.symbols;
        symbol_freq_blocks[b].resize(local_model.symbols.size());
        for (size_t i = 0; i < local_model.symbols.size(); ++i)
            symbol_freq_blocks[b][i] = local_model.symbols_freq[i];
        total_preterm_blocks[b] = local_model.total_preterm;
    }

    // 3. 合并所有块的统计数据到主模型（流水线归并）
    // 只需串行遍历每个块，合并到全局哈希表
    unordered_map<string, int> pt_hash;
    unordered_map<string, int> letter_hash, digit_hash, symbol_hash;
    auto pt2str = [](const PT& pt) {
        string s;
        for (const auto& seg : pt.content) {
            s += to_string(seg.type) + "-" + to_string(seg.length) + "|";
        }
        return s;
    };
    auto seg2str = [](const segment& seg) {
        return to_string(seg.type) + "-" + to_string(seg.length);
    };

    for (size_t b = 0; b < num_blocks; ++b) {
        // PT
        for (size_t i = 0; i < pt_blocks[b].size(); ++i) {
            string key = pt2str(pt_blocks[b][i]);
            auto it = pt_hash.find(key);
            if (it == pt_hash.end()) {
                int new_id = GetNextPretermID();
                preterminals.push_back(pt_blocks[b][i]);
                preterm_freq[new_id] = pt_freq_blocks[b][i];
                pt_hash[key] = new_id;
            } else {
                preterm_freq[it->second] += pt_freq_blocks[b][i];
            }
        }
        // letters
        for (size_t i = 0; i < letter_blocks[b].size(); ++i) {
            string key = seg2str(letter_blocks[b][i]);
            auto it = letter_hash.find(key);
            if (it == letter_hash.end()) {
                int new_id = GetNextLettersID();
                letters.push_back(letter_blocks[b][i]);
                letters_freq[new_id] = letter_freq_blocks[b][i];
                letter_hash[key] = new_id;
            } else {
                int id = it->second;
                letters_freq[id] += letter_freq_blocks[b][i];
                for (const auto& kv : letter_blocks[b][i].values) {
                    const string& val = kv.first;
                    int idx = kv.second;
                    if (letters[id].values.find(val) == letters[id].values.end()) {
                        int new_idx = letters[id].values.size();
                        letters[id].values[val] = new_idx;
                        letters[id].freqs[new_idx] = letter_blocks[b][i].freqs.at(idx);
                    } else {
                        int exist_idx = letters[id].values[val];
                        letters[id].freqs[exist_idx] += letter_blocks[b][i].freqs.at(idx);
                    }
                }
            }
        }
        // digits
        for (size_t i = 0; i < digit_blocks[b].size(); ++i) {
            string key = seg2str(digit_blocks[b][i]);
            auto it = digit_hash.find(key);
            if (it == digit_hash.end()) {
                int new_id = GetNextDigitsID();
                digits.push_back(digit_blocks[b][i]);
                digits_freq[new_id] = digit_freq_blocks[b][i];
                digit_hash[key] = new_id;
            } else {
                int id = it->second;
                digits_freq[id] += digit_freq_blocks[b][i];
                for (const auto& kv : digit_blocks[b][i].values) {
                    const string& val = kv.first;
                    int idx = kv.second;
                    if (digits[id].values.find(val) == digits[id].values.end()) {
                        int new_idx = digits[id].values.size();
                        digits[id].values[val] = new_idx;
                        digits[id].freqs[new_idx] = digit_blocks[b][i].freqs.at(idx);
                    } else {
                        int exist_idx = digits[id].values[val];
                        digits[id].freqs[exist_idx] += digit_blocks[b][i].freqs.at(idx);
                    }
                }
            }
        }
        // symbols
        for (size_t i = 0; i < symbol_blocks[b].size(); ++i) {
            string key = seg2str(symbol_blocks[b][i]);
            auto it = symbol_hash.find(key);
            if (it == symbol_hash.end()) {
                int new_id = GetNextSymbolsID();
                symbols.push_back(symbol_blocks[b][i]);
                symbols_freq[new_id] = symbol_freq_blocks[b][i];
                symbol_hash[key] = new_id;
            } else {
                int id = it->second;
                symbols_freq[id] += symbol_freq_blocks[b][i];
                for (const auto& kv : symbol_blocks[b][i].values) {
                    const string& val = kv.first;
                    int idx = kv.second;
                    if (symbols[id].values.find(val) == symbols[id].values.end()) {
                        int new_idx = symbols[id].values.size();
                        symbols[id].values[val] = new_idx;
                        symbols[id].freqs[new_idx] = symbol_blocks[b][i].freqs.at(idx);
                    } else {
                        int exist_idx = symbols[id].values[val];
                        symbols[id].freqs[exist_idx] += symbol_blocks[b][i].freqs.at(idx);
                    }
                }
            }
        }
    }
    total_preterm = 0;
    for (size_t b = 0; b < num_blocks; ++b) total_preterm += total_preterm_blocks[b];

    cout << "Training phase 1 complete. Total passwords processed: " << lines << endl;
}

/// @brief 在模型中找到一个PT的统计数据
/// @param pt 需要查找的PT
/// @return 目标PT在模型中的对应下标
int model::FindPT(PT pt)
{
    for (int id = 0; id < preterminals.size(); id += 1)
    {
        if (preterminals[id].content.size() != pt.content.size())
        {
            continue;
        }
        else
        {
            bool equal_flag = true;
            for (int idx = 0; idx < preterminals[id].content.size(); idx += 1)
            {
                if (preterminals[id].content[idx].type != pt.content[idx].type || preterminals[id].content[idx].length != pt.content[idx].length)
                {
                    equal_flag = false;
                    break;
                }
            }
            if (equal_flag == true)
            {
                return id;
            }
        }
    }
    return -1;
}

/// @brief 在模型中找到一个letter segment的统计数据
/// @param seg 要找的letter segment
/// @return 目标letter segment的对应下标
int model::FindLetter(segment seg)
{
    for (int id = 0; id < letters.size(); id += 1)
    {
        if (letters[id].length == seg.length)
        {
            return id;
        }
    }
    return -1;
}

/// @brief 在模型中找到一个digit segment的统计数据
/// @param seg 要找的digit segment
/// @return 目标digit segment的对应下标
int model::FindDigit(segment seg)
{
    for (int id = 0; id < digits.size(); id += 1)
    {
        if (digits[id].length == seg.length)
        {
            return id;
        }
    }
    return -1;
}

int model::FindSymbol(segment seg)
{
    for (int id = 0; id < symbols.size(); id += 1)
    {
        if (symbols[id].length == seg.length)
        {
            return id;
        }
    }
    return -1;
}

void PT::insert(segment seg)
{
    content.emplace_back(seg);
}

void segment::insert(string value)
{
    if (values.find(value) == values.end())
    {
        values[value] = values.size();
        freqs[values[value]] = 1;
    }
    else
    {
        freqs[values[value]] += 1;
    }
}


void segment::order()
{
    ordered_values.reserve(values.size());
    for (const auto& value : values)
    {
        ordered_values.emplace_back(value.first);
    }
    sort(ordered_values.begin(), ordered_values.end(),
              [this](const string &a, const string &b)
              {
                  return freqs.at(values.at(a)) > freqs.at(values.at(b));
              });

    ordered_freqs.reserve(ordered_values.size());
    total_freq = 0;
    for (const string &val : ordered_values)
    {
        ordered_freqs.emplace_back(freqs.at(values.at(val)));
        total_freq += freqs.at(values.at(val));
    }
}

void model::parse(string pw)
{
    PT pt;
    string curr_part = "";
    int curr_type = 0; // 0: 未设置, 1: 字母, 2: 数字, 3: 特殊字符
    // 请学会使用这种方式写for循环：for (auto it : iterable)
    // 相信我，以后你会用上的。You're welcome :)
    for (char ch : pw)
    {
        if (isalpha(ch))
        {
            if (curr_type != 1)
            {
                if (curr_type == 2)
                {
                    segment seg(curr_type, curr_part.length());
                    if (FindDigit(seg) == -1)
                    {
                        int id = GetNextDigitsID();
                        digits.emplace_back(seg);
                        digits[id].insert(curr_part);
                        digits_freq[id] = 1;
                    }
                    else
                    {
                        int id = FindDigit(seg);
                        digits_freq[id] += 1;
                        digits[id].insert(curr_part);
                    }
                    curr_part.clear();
                    pt.insert(seg);
                }
                else if (curr_type == 3)
                {
                    segment seg(curr_type, curr_part.length());
                    if (FindSymbol(seg) == -1)
                    {
                        int id = GetNextSymbolsID();
                        symbols.emplace_back(seg);
                        symbols_freq[id] = 1;
                        symbols[id].insert(curr_part);
                    }
                    else
                    {
                        int id = FindSymbol(seg);
                        symbols_freq[id] += 1;
                        symbols[id].insert(curr_part);
                    }
                    curr_part.clear();
                    pt.insert(seg);
                }
            }
            curr_type = 1;
            curr_part += ch;
        }
        else if (isdigit(ch))
        {
            if (curr_type != 2)
            {
                if (curr_type == 1)
                {
                    segment seg(curr_type, curr_part.length());
                    if (FindLetter(seg) == -1)
                    {
                        int id = GetNextLettersID();
                        letters.emplace_back(seg);
                        letters_freq[id] = 1;
                        letters[id].insert(curr_part);
                    }
                    else
                    {
                        int id = FindLetter(seg);
                        letters_freq[id] += 1;
                        letters[id].insert(curr_part);
                    }
                    curr_part.clear();
                    pt.insert(seg);
                }
                else if (curr_type == 3)
                {
                    segment seg(curr_type, curr_part.length());
                    if (FindSymbol(seg) == -1)
                    {
                        int id = GetNextSymbolsID();
                        symbols.emplace_back(seg);
                        symbols_freq[id] = 1;
                        symbols[id].insert(curr_part);
                    }
                    else
                    {
                        int id = FindSymbol(seg);
                        symbols_freq[id] += 1;
                        symbols[id].insert(curr_part);
                    }
                    curr_part.clear();
                    pt.insert(seg);
                }
            }
            curr_type = 2;
            curr_part += ch;
        }
        else
        {
            if (curr_type != 3)
            {
                if (curr_type == 1)
                {
                    segment seg(curr_type, curr_part.length());
                    if (FindLetter(seg) == -1)
                    {
                        int id = GetNextLettersID();
                        letters.emplace_back(seg);
                        letters_freq[id] = 1;
                        letters[id].insert(curr_part);
                    }
                    else
                    {
                        int id = FindLetter(seg);
                        letters_freq[id] += 1;
                        letters[id].insert(curr_part);
                    }
                    curr_part.clear();
                    pt.insert(seg);
                }
                else if (curr_type == 2)
                {
                    segment seg(curr_type, curr_part.length());
                    if (FindDigit(seg) == -1)
                    {
                        int id = GetNextDigitsID();
                        digits.emplace_back(seg);
                        digits_freq[id] = 1;
                        digits[id].insert(curr_part);
                    }
                    else
                    {
                        int id = FindDigit(seg);
                        digits_freq[id] += 1;
                        digits[id].insert(curr_part);
                    }
                    curr_part.clear();
                    pt.insert(seg);
                }
            }
            curr_type = 3;
            curr_part += ch;
        }
    }
    if (!curr_part.empty())
    {
        if (curr_type == 1)
        {
            segment seg(curr_type, curr_part.length());
            if (FindLetter(seg) == -1)
            {
                int id = GetNextLettersID();
                letters.emplace_back(seg);
                letters_freq[id] = 1;
                letters[id].insert(curr_part);
            }
            else
            {
                int id = FindLetter(seg);
                letters_freq[id] += 1;
                letters[id].insert(curr_part);
            }
            curr_part.clear();
            pt.insert(seg);
        }
        else if (curr_type == 2)
        {
            segment seg(curr_type, curr_part.length());
            if (FindDigit(seg) == -1)
            {
                int id = GetNextDigitsID();
                digits.emplace_back(seg);
                digits_freq[id] = 1;
                digits[id].insert(curr_part);
            }
            else
            {
                int id = FindDigit(seg);
                digits_freq[id] += 1;
                digits[id].insert(curr_part);
            }
            curr_part.clear();
            pt.insert(seg);
        }
        else
        {
            segment seg(curr_type, curr_part.length());
            if (FindSymbol(seg) == -1)
            {
                int id = GetNextSymbolsID();
                symbols.emplace_back(seg);
                symbols_freq[id] = 1;
                symbols[id].insert(curr_part);
            }
            else
            {
                int id = FindSymbol(seg);
                symbols_freq[id] += 1;
                symbols[id].insert(curr_part);
            }
            curr_part.clear();
            pt.insert(seg);
        }
    }
    // pt.PrintPT();
    // cout<<endl;
    // cout << FindPT(pt) << endl;
    total_preterm += 1;
    if (FindPT(pt) == -1)
    {
        for (int i = 0; i < pt.content.size(); i += 1)
        {
            pt.curr_indices.emplace_back(0);
        }
        int id = GetNextPretermID();
        // cout << id << endl;
        preterminals.emplace_back(pt);
        preterm_freq[id] = 1;
    }
    else
    {
        int id = FindPT(pt);
        // cout << id << endl;
        preterm_freq[id] += 1;
    }
}

void segment::PrintSeg()
{
    if (type == 1)
    {
        cout << "L" << length;
    }
    if (type == 2)
    {
        cout << "D" << length;
    }
    if (type == 3)
    {
        cout << "S" << length;
    }
}

void segment::PrintValues()
{
    // order();
    for (string iter : ordered_values)
    {
        cout << iter << " freq:" << freqs[values[iter]] << endl;
    }
}

void PT::PrintPT()
{
    for (auto iter : content)
    {
        iter.PrintSeg();
    }
}

void model::print()
{
    cout << "preterminals:" << endl;
    for (int i = 0; i < preterminals.size(); i += 1)
    {
        preterminals[i].PrintPT();
        // cout << preterminals[i].curr_indices.size() << endl;
        cout << " freq:" << preterm_freq[i];
        cout << endl;
    }
    // order();
    for (auto iter : ordered_pts)
    {
        iter.PrintPT();
        cout << " freq:" << preterm_freq[FindPT(iter)];
        cout << endl;
    }
    cout << "segments:" << endl;
    for (int i = 0; i < letters.size(); i += 1)
    {
        letters[i].PrintSeg();
        // letters[i].PrintValues();
        cout << " freq:" << letters_freq[i];
        cout << endl;
    }
    for (int i = 0; i < digits.size(); i += 1)
    {
        digits[i].PrintSeg();
        // digits[i].PrintValues();
        cout << " freq:" << digits_freq[i];
        cout << endl;
    }
    for (int i = 0; i < symbols.size(); i += 1)
    {
        symbols[i].PrintSeg();
        // symbols[i].PrintValues();
        cout << " freq:" << symbols_freq[i];
        cout << endl;
    }
}

bool compareByPretermProb(const PT& a, const PT& b) {
    return a.preterm_prob > b.preterm_prob;  // 降序排序
}

void model::order()
{
    cout << "Training phase 2: Ordering segment values and PTs..." << endl;

    // 1. 并行填充ordered_pts和计算preterm_prob
    ordered_pts.resize(preterminals.size());
    #pragma omp parallel for
    for (int i = 0; i < (int)preterminals.size(); ++i) {
        PT pt = preterminals[i];
        pt.preterm_prob = float(preterm_freq[FindPT(pt)]) / total_preterm;
        ordered_pts[i] = pt;
    }

    cout << "total pts" << ordered_pts.size() << endl;
    // 排序本身一般不需要并行，sort已很快
    sort(ordered_pts.begin(), ordered_pts.end(), compareByPretermProb);

    // 2. 并行对letters/digits/symbols排序
    cout << "Ordering letters" << endl;
    #pragma omp parallel for
    for (int i = 0; i < (int)letters.size(); ++i) {
        letters[i].order();
    }
    cout << "Ordering digits" << endl;
    #pragma omp parallel for
    for (int i = 0; i < (int)digits.size(); ++i) {
        digits[i].order();
    }
    cout << "ordering symbols" << endl;
    #pragma omp parallel for
    for (int i = 0; i < (int)symbols.size(); ++i) {
        symbols[i].order();
    }
    cout << "Training phase 2 complete." << endl;
=======
#include "PCFG.h"
#include <fstream>
#include <cctype>
#include <algorithm>
#include <omp.h>
#include <vector>
#include <unordered_map>

// 这个文件里面的各函数你都不需要完全理解，甚至根本不需要看
// 从学术价值上讲，加速模型的训练过程是一个没什么价值的问题，因为我们一般假定统计学模型的训练成本较低
// 但是，假如你是一个投稿时顶着ddl做实验的倒霉研究生/实习生，提高训练速度就可以大幅节省你的时间了
// 所以如果你愿意，也可以尝试用多线程加速训练过程

/**
 * 怎么加速PCFG训练过程？据助教所知，没有公开文献提出过有效的加速方法（因为这么做基本无学术价值）
 * 
 * 但是统计学模型好就好在其数据是可加的。例如，假如我把数据集拆分成4个部分，并行训练4个不同的模型。
 * 然后我可以直接将四个模型的统计数据进行简单加和，就得到了和串行训练相同的模型了。
 * 
 * 说起来容易，做起来不一定容易，你可能会碰到一系列具体的工程问题。如果你决定加速训练过程，祝你好运！
 * 
 */

// 训练的wrapper，实际上就是读取训练集
void model::train(string path)
{
    string pw;
    ifstream train_set(path);
    int lines = 0;
    cout<<"Training..."<<endl;
    cout<<"Training phase 1: reading and parsing passwords..."<<endl;

    // 1. 先把所有密码读入内存
    vector<string> passwords;
    passwords.reserve(3000000);
    while (train_set >> pw)
    {
        lines += 1;
        if (lines % 100000 == 0)
        {
            cout <<"Lines processed: "<< lines << endl;
            if (lines > 3000000)
            {
                break;
            }
        }
        passwords.push_back(pw);
    }

    int n_threads = omp_get_max_threads();
    vector<model> local_models(n_threads);

    // 2. 流水线并行处理：分块读取、分块解析、分块合并
    const size_t block_size = 100000; // 可根据内存和核数调整
    size_t total = passwords.size();
    size_t num_blocks = (total + block_size - 1) / block_size;

    // 每个线程有自己的本地桶，减少同步
    vector<vector<PT>> pt_blocks(num_blocks);
    vector<vector<int>> pt_freq_blocks(num_blocks);
    vector<vector<segment>> letter_blocks(num_blocks), digit_blocks(num_blocks), symbol_blocks(num_blocks);
    vector<vector<int>> letter_freq_blocks(num_blocks), digit_freq_blocks(num_blocks), symbol_freq_blocks(num_blocks);
    vector<int> total_preterm_blocks(num_blocks, 0);

    #pragma omp parallel for schedule(dynamic)
    for (size_t b = 0; b < num_blocks; ++b) {
        size_t start = b * block_size;
        size_t end = min(start + block_size, total);
        model local_model;
        local_model.preterminals.reserve(block_size / 2);
        local_model.letters.reserve(128);
        local_model.digits.reserve(128);
        local_model.symbols.reserve(128);
        for (size_t i = start; i < end; ++i) {
            local_model.parse(passwords[i]);
        }
        // 收集本块数据
        pt_blocks[b] = local_model.preterminals;
        pt_freq_blocks[b].resize(local_model.preterminals.size());
        for (size_t i = 0; i < local_model.preterminals.size(); ++i)
            pt_freq_blocks[b][i] = local_model.preterm_freq[i];
        letter_blocks[b] = local_model.letters;
        letter_freq_blocks[b].resize(local_model.letters.size());
        for (size_t i = 0; i < local_model.letters.size(); ++i)
            letter_freq_blocks[b][i] = local_model.letters_freq[i];
        digit_blocks[b] = local_model.digits;
        digit_freq_blocks[b].resize(local_model.digits.size());
        for (size_t i = 0; i < local_model.digits.size(); ++i)
            digit_freq_blocks[b][i] = local_model.digits_freq[i];
        symbol_blocks[b] = local_model.symbols;
        symbol_freq_blocks[b].resize(local_model.symbols.size());
        for (size_t i = 0; i < local_model.symbols.size(); ++i)
            symbol_freq_blocks[b][i] = local_model.symbols_freq[i];
        total_preterm_blocks[b] = local_model.total_preterm;
    }

    // 3. 合并所有块的统计数据到主模型（流水线归并）
    // 只需串行遍历每个块，合并到全局哈希表
    unordered_map<string, int> pt_hash;
    unordered_map<string, int> letter_hash, digit_hash, symbol_hash;
    auto pt2str = [](const PT& pt) {
        string s;
        for (const auto& seg : pt.content) {
            s += to_string(seg.type) + "-" + to_string(seg.length) + "|";
        }
        return s;
    };
    auto seg2str = [](const segment& seg) {
        return to_string(seg.type) + "-" + to_string(seg.length);
    };

    for (size_t b = 0; b < num_blocks; ++b) {
        // PT
        for (size_t i = 0; i < pt_blocks[b].size(); ++i) {
            string key = pt2str(pt_blocks[b][i]);
            auto it = pt_hash.find(key);
            if (it == pt_hash.end()) {
                int new_id = GetNextPretermID();
                preterminals.push_back(pt_blocks[b][i]);
                preterm_freq[new_id] = pt_freq_blocks[b][i];
                pt_hash[key] = new_id;
            } else {
                preterm_freq[it->second] += pt_freq_blocks[b][i];
            }
        }
        // letters
        for (size_t i = 0; i < letter_blocks[b].size(); ++i) {
            string key = seg2str(letter_blocks[b][i]);
            auto it = letter_hash.find(key);
            if (it == letter_hash.end()) {
                int new_id = GetNextLettersID();
                letters.push_back(letter_blocks[b][i]);
                letters_freq[new_id] = letter_freq_blocks[b][i];
                letter_hash[key] = new_id;
            } else {
                int id = it->second;
                letters_freq[id] += letter_freq_blocks[b][i];
                for (const auto& kv : letter_blocks[b][i].values) {
                    const string& val = kv.first;
                    int idx = kv.second;
                    if (letters[id].values.find(val) == letters[id].values.end()) {
                        int new_idx = letters[id].values.size();
                        letters[id].values[val] = new_idx;
                        letters[id].freqs[new_idx] = letter_blocks[b][i].freqs.at(idx);
                    } else {
                        int exist_idx = letters[id].values[val];
                        letters[id].freqs[exist_idx] += letter_blocks[b][i].freqs.at(idx);
                    }
                }
            }
        }
        // digits
        for (size_t i = 0; i < digit_blocks[b].size(); ++i) {
            string key = seg2str(digit_blocks[b][i]);
            auto it = digit_hash.find(key);
            if (it == digit_hash.end()) {
                int new_id = GetNextDigitsID();
                digits.push_back(digit_blocks[b][i]);
                digits_freq[new_id] = digit_freq_blocks[b][i];
                digit_hash[key] = new_id;
            } else {
                int id = it->second;
                digits_freq[id] += digit_freq_blocks[b][i];
                for (const auto& kv : digit_blocks[b][i].values) {
                    const string& val = kv.first;
                    int idx = kv.second;
                    if (digits[id].values.find(val) == digits[id].values.end()) {
                        int new_idx = digits[id].values.size();
                        digits[id].values[val] = new_idx;
                        digits[id].freqs[new_idx] = digit_blocks[b][i].freqs.at(idx);
                    } else {
                        int exist_idx = digits[id].values[val];
                        digits[id].freqs[exist_idx] += digit_blocks[b][i].freqs.at(idx);
                    }
                }
            }
        }
        // symbols
        for (size_t i = 0; i < symbol_blocks[b].size(); ++i) {
            string key = seg2str(symbol_blocks[b][i]);
            auto it = symbol_hash.find(key);
            if (it == symbol_hash.end()) {
                int new_id = GetNextSymbolsID();
                symbols.push_back(symbol_blocks[b][i]);
                symbols_freq[new_id] = symbol_freq_blocks[b][i];
                symbol_hash[key] = new_id;
            } else {
                int id = it->second;
                symbols_freq[id] += symbol_freq_blocks[b][i];
                for (const auto& kv : symbol_blocks[b][i].values) {
                    const string& val = kv.first;
                    int idx = kv.second;
                    if (symbols[id].values.find(val) == symbols[id].values.end()) {
                        int new_idx = symbols[id].values.size();
                        symbols[id].values[val] = new_idx;
                        symbols[id].freqs[new_idx] = symbol_blocks[b][i].freqs.at(idx);
                    } else {
                        int exist_idx = symbols[id].values[val];
                        symbols[id].freqs[exist_idx] += symbol_blocks[b][i].freqs.at(idx);
                    }
                }
            }
        }
    }
    total_preterm = 0;
    for (size_t b = 0; b < num_blocks; ++b) total_preterm += total_preterm_blocks[b];

    cout << "Training phase 1 complete. Total passwords processed: " << lines << endl;
}

/// @brief 在模型中找到一个PT的统计数据
/// @param pt 需要查找的PT
/// @return 目标PT在模型中的对应下标
int model::FindPT(PT pt)
{
    for (int id = 0; id < preterminals.size(); id += 1)
    {
        if (preterminals[id].content.size() != pt.content.size())
        {
            continue;
        }
        else
        {
            bool equal_flag = true;
            for (int idx = 0; idx < preterminals[id].content.size(); idx += 1)
            {
                if (preterminals[id].content[idx].type != pt.content[idx].type || preterminals[id].content[idx].length != pt.content[idx].length)
                {
                    equal_flag = false;
                    break;
                }
            }
            if (equal_flag == true)
            {
                return id;
            }
        }
    }
    return -1;
}

/// @brief 在模型中找到一个letter segment的统计数据
/// @param seg 要找的letter segment
/// @return 目标letter segment的对应下标
int model::FindLetter(segment seg)
{
    for (int id = 0; id < letters.size(); id += 1)
    {
        if (letters[id].length == seg.length)
        {
            return id;
        }
    }
    return -1;
}

/// @brief 在模型中找到一个digit segment的统计数据
/// @param seg 要找的digit segment
/// @return 目标digit segment的对应下标
int model::FindDigit(segment seg)
{
    for (int id = 0; id < digits.size(); id += 1)
    {
        if (digits[id].length == seg.length)
        {
            return id;
        }
    }
    return -1;
}

int model::FindSymbol(segment seg)
{
    for (int id = 0; id < symbols.size(); id += 1)
    {
        if (symbols[id].length == seg.length)
        {
            return id;
        }
    }
    return -1;
}

void PT::insert(segment seg)
{
    content.emplace_back(seg);
}

void segment::insert(string value)
{
    if (values.find(value) == values.end())
    {
        values[value] = values.size();
        freqs[values[value]] = 1;
    }
    else
    {
        freqs[values[value]] += 1;
    }
}


void segment::order()
{
    ordered_values.reserve(values.size());
    for (const auto& value : values)
    {
        ordered_values.emplace_back(value.first);
    }
    sort(ordered_values.begin(), ordered_values.end(),
              [this](const string &a, const string &b)
              {
                  return freqs.at(values.at(a)) > freqs.at(values.at(b));
              });

    ordered_freqs.reserve(ordered_values.size());
    total_freq = 0;
    for (const string &val : ordered_values)
    {
        ordered_freqs.emplace_back(freqs.at(values.at(val)));
        total_freq += freqs.at(values.at(val));
    }
}

void model::parse(string pw)
{
    PT pt;
    string curr_part = "";
    int curr_type = 0; // 0: 未设置, 1: 字母, 2: 数字, 3: 特殊字符
    // 请学会使用这种方式写for循环：for (auto it : iterable)
    // 相信我，以后你会用上的。You're welcome :)
    for (char ch : pw)
    {
        if (isalpha(ch))
        {
            if (curr_type != 1)
            {
                if (curr_type == 2)
                {
                    segment seg(curr_type, curr_part.length());
                    if (FindDigit(seg) == -1)
                    {
                        int id = GetNextDigitsID();
                        digits.emplace_back(seg);
                        digits[id].insert(curr_part);
                        digits_freq[id] = 1;
                    }
                    else
                    {
                        int id = FindDigit(seg);
                        digits_freq[id] += 1;
                        digits[id].insert(curr_part);
                    }
                    curr_part.clear();
                    pt.insert(seg);
                }
                else if (curr_type == 3)
                {
                    segment seg(curr_type, curr_part.length());
                    if (FindSymbol(seg) == -1)
                    {
                        int id = GetNextSymbolsID();
                        symbols.emplace_back(seg);
                        symbols_freq[id] = 1;
                        symbols[id].insert(curr_part);
                    }
                    else
                    {
                        int id = FindSymbol(seg);
                        symbols_freq[id] += 1;
                        symbols[id].insert(curr_part);
                    }
                    curr_part.clear();
                    pt.insert(seg);
                }
            }
            curr_type = 1;
            curr_part += ch;
        }
        else if (isdigit(ch))
        {
            if (curr_type != 2)
            {
                if (curr_type == 1)
                {
                    segment seg(curr_type, curr_part.length());
                    if (FindLetter(seg) == -1)
                    {
                        int id = GetNextLettersID();
                        letters.emplace_back(seg);
                        letters_freq[id] = 1;
                        letters[id].insert(curr_part);
                    }
                    else
                    {
                        int id = FindLetter(seg);
                        letters_freq[id] += 1;
                        letters[id].insert(curr_part);
                    }
                    curr_part.clear();
                    pt.insert(seg);
                }
                else if (curr_type == 3)
                {
                    segment seg(curr_type, curr_part.length());
                    if (FindSymbol(seg) == -1)
                    {
                        int id = GetNextSymbolsID();
                        symbols.emplace_back(seg);
                        symbols_freq[id] = 1;
                        symbols[id].insert(curr_part);
                    }
                    else
                    {
                        int id = FindSymbol(seg);
                        symbols_freq[id] += 1;
                        symbols[id].insert(curr_part);
                    }
                    curr_part.clear();
                    pt.insert(seg);
                }
            }
            curr_type = 2;
            curr_part += ch;
        }
        else
        {
            if (curr_type != 3)
            {
                if (curr_type == 1)
                {
                    segment seg(curr_type, curr_part.length());
                    if (FindLetter(seg) == -1)
                    {
                        int id = GetNextLettersID();
                        letters.emplace_back(seg);
                        letters_freq[id] = 1;
                        letters[id].insert(curr_part);
                    }
                    else
                    {
                        int id = FindLetter(seg);
                        letters_freq[id] += 1;
                        letters[id].insert(curr_part);
                    }
                    curr_part.clear();
                    pt.insert(seg);
                }
                else if (curr_type == 2)
                {
                    segment seg(curr_type, curr_part.length());
                    if (FindDigit(seg) == -1)
                    {
                        int id = GetNextDigitsID();
                        digits.emplace_back(seg);
                        digits_freq[id] = 1;
                        digits[id].insert(curr_part);
                    }
                    else
                    {
                        int id = FindDigit(seg);
                        digits_freq[id] += 1;
                        digits[id].insert(curr_part);
                    }
                    curr_part.clear();
                    pt.insert(seg);
                }
            }
            curr_type = 3;
            curr_part += ch;
        }
    }
    if (!curr_part.empty())
    {
        if (curr_type == 1)
        {
            segment seg(curr_type, curr_part.length());
            if (FindLetter(seg) == -1)
            {
                int id = GetNextLettersID();
                letters.emplace_back(seg);
                letters_freq[id] = 1;
                letters[id].insert(curr_part);
            }
            else
            {
                int id = FindLetter(seg);
                letters_freq[id] += 1;
                letters[id].insert(curr_part);
            }
            curr_part.clear();
            pt.insert(seg);
        }
        else if (curr_type == 2)
        {
            segment seg(curr_type, curr_part.length());
            if (FindDigit(seg) == -1)
            {
                int id = GetNextDigitsID();
                digits.emplace_back(seg);
                digits_freq[id] = 1;
                digits[id].insert(curr_part);
            }
            else
            {
                int id = FindDigit(seg);
                digits_freq[id] += 1;
                digits[id].insert(curr_part);
            }
            curr_part.clear();
            pt.insert(seg);
        }
        else
        {
            segment seg(curr_type, curr_part.length());
            if (FindSymbol(seg) == -1)
            {
                int id = GetNextSymbolsID();
                symbols.emplace_back(seg);
                symbols_freq[id] = 1;
                symbols[id].insert(curr_part);
            }
            else
            {
                int id = FindSymbol(seg);
                symbols_freq[id] += 1;
                symbols[id].insert(curr_part);
            }
            curr_part.clear();
            pt.insert(seg);
        }
    }
    // pt.PrintPT();
    // cout<<endl;
    // cout << FindPT(pt) << endl;
    total_preterm += 1;
    if (FindPT(pt) == -1)
    {
        for (int i = 0; i < pt.content.size(); i += 1)
        {
            pt.curr_indices.emplace_back(0);
        }
        int id = GetNextPretermID();
        // cout << id << endl;
        preterminals.emplace_back(pt);
        preterm_freq[id] = 1;
    }
    else
    {
        int id = FindPT(pt);
        // cout << id << endl;
        preterm_freq[id] += 1;
    }
}

void segment::PrintSeg()
{
    if (type == 1)
    {
        cout << "L" << length;
    }
    if (type == 2)
    {
        cout << "D" << length;
    }
    if (type == 3)
    {
        cout << "S" << length;
    }
}

void segment::PrintValues()
{
    // order();
    for (string iter : ordered_values)
    {
        cout << iter << " freq:" << freqs[values[iter]] << endl;
    }
}

void PT::PrintPT()
{
    for (auto iter : content)
    {
        iter.PrintSeg();
    }
}

void model::print()
{
    cout << "preterminals:" << endl;
    for (int i = 0; i < preterminals.size(); i += 1)
    {
        preterminals[i].PrintPT();
        // cout << preterminals[i].curr_indices.size() << endl;
        cout << " freq:" << preterm_freq[i];
        cout << endl;
    }
    // order();
    for (auto iter : ordered_pts)
    {
        iter.PrintPT();
        cout << " freq:" << preterm_freq[FindPT(iter)];
        cout << endl;
    }
    cout << "segments:" << endl;
    for (int i = 0; i < letters.size(); i += 1)
    {
        letters[i].PrintSeg();
        // letters[i].PrintValues();
        cout << " freq:" << letters_freq[i];
        cout << endl;
    }
    for (int i = 0; i < digits.size(); i += 1)
    {
        digits[i].PrintSeg();
        // digits[i].PrintValues();
        cout << " freq:" << digits_freq[i];
        cout << endl;
    }
    for (int i = 0; i < symbols.size(); i += 1)
    {
        symbols[i].PrintSeg();
        // symbols[i].PrintValues();
        cout << " freq:" << symbols_freq[i];
        cout << endl;
    }
}

bool compareByPretermProb(const PT& a, const PT& b) {
    return a.preterm_prob > b.preterm_prob;  // 降序排序
}

void model::order()
{
    cout << "Training phase 2: Ordering segment values and PTs..." << endl;

    // 1. 并行填充ordered_pts和计算preterm_prob
    ordered_pts.resize(preterminals.size());
    #pragma omp parallel for
    for (int i = 0; i < (int)preterminals.size(); ++i) {
        PT pt = preterminals[i];
        pt.preterm_prob = float(preterm_freq[FindPT(pt)]) / total_preterm;
        ordered_pts[i] = pt;
    }

    cout << "total pts" << ordered_pts.size() << endl;
    // 排序本身一般不需要并行，sort已很快
    sort(ordered_pts.begin(), ordered_pts.end(), compareByPretermProb);

    // 2. 并行对letters/digits/symbols排序
    cout << "Ordering letters" << endl;
    #pragma omp parallel for
    for (int i = 0; i < (int)letters.size(); ++i) {
        letters[i].order();
    }
    cout << "Ordering digits" << endl;
    #pragma omp parallel for
    for (int i = 0; i < (int)digits.size(); ++i) {
        digits[i].order();
    }
    cout << "ordering symbols" << endl;
    #pragma omp parallel for
    for (int i = 0; i < (int)symbols.size(); ++i) {
        symbols[i].order();
    }
    cout << "Training phase 2 complete." << endl;
>>>>>>> 414c772c5f5b372157d1cd655d0951d27e7baec0
}