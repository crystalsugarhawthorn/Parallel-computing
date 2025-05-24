#include <iostream>
#include <string>
#include <cstring>
#include <arm_neon.h>

using namespace std;

// 定义了Byte，便于使用
typedef unsigned char Byte;
// 定义了32比特
typedef unsigned int bit32;

// MD5的一系列参数。参数是固定的，其实你不需要看懂这些
#define s11 7
#define s12 12
#define s13 17
#define s14 22
#define s21 5
#define s22 9
#define s23 14
#define s24 20
#define s31 4
#define s32 11
#define s33 16
#define s34 23
#define s41 6
#define s42 10
#define s43 15
#define s44 21

/**
 * @Basic MD5 functions.
 *
 * @param there bit32.
 *
 * @return one bit32.
 */
// 定义了一系列MD5中的具体函数
// 这四个计算函数是需要你进行SIMD并行化的
// 可以看到，FGHI四个函数都涉及一系列位运算，在数据上是对齐的，非常容易实现SIMD的并行化

//--------串行部分-------//
#define F(x, y, z) (((x) & (y)) | ((~x) & (z)))
#define G(x, y, z) (((x) & (z)) | ((y) & (~z)))
#define H(x, y, z) ((x) ^ (y) ^ (z))
#define I(x, y, z) ((y) ^ ((x) | (~z)))


//-------SIMD并行化部分-------//
/*
inline uint32x4_t F_neon(uint32x4_t x, uint32x4_t y, uint32x4_t z) {
  uint32x4_t term1 = vandq_u32(x, y);        // x & y
  uint32x4_t term2 = vandq_u32(vmvnq_u32(x), z); // ~x & z
  return vorrq_u32(term1, term2);            // 按位或
}
inline uint32x4_t G_neon(uint32x4_t x, uint32x4_t y, uint32x4_t z) {
  uint32x4_t term1 = vandq_u32(x, z);        // x & z
  uint32x4_t term2 = vandq_u32(y, vmvnq_u32(z)); // y & ~z
  return vorrq_u32(term1, term2);            // 按位或
}
inline uint32x4_t H_neon(uint32x4_t x, uint32x4_t y, uint32x4_t z) {
  uint32x4_t tmp = veorq_u32(x, y);          // x ^ y
  return veorq_u32(tmp, z);                  // 结果再与z异或
}
inline uint32x4_t I_neon(uint32x4_t x, uint32x4_t y, uint32x4_t z) {
  uint32x4_t not_z = vmvnq_u32(z);           // ~z
  uint32x4_t or_term = vorrq_u32(x, not_z);  // x | ~z
  return veorq_u32(y, or_term);              // y异或结果
}
*/

// 定义宏（直接展开为NEON指令组合）
// F(x,y,z) = (x & y) | (~x & z)
#define F_NEON(x, y, z) \
    vorrq_u32(vandq_u32((x), (y)), vandq_u32(vmvnq_u32((x)), (z)))

// G(x,y,z) = (x & z) | (y & ~z)
#define G_NEON(x, y, z) \
    vorrq_u32(vandq_u32((x), (z)), vandq_u32((y), vmvnq_u32((z))))

// H(x,y,z) = x ^ y ^ z
#define H_NEON(x, y, z) \
    veorq_u32(veorq_u32((x), (y)), (z))

// I(x,y,z) = y ^ (x | ~z)
#define I_NEON(x, y, z) \
    veorq_u32((y), vorrq_u32((x), vmvnq_u32((z))))


/**
 * @Rotate Left.
 *
 * @param {num} the raw number.
 *
 * @param {n} rotate left n.
 *
 * @return the number after rotated left.
 */
// 定义了一系列MD5中的具体函数
// 这五个计算函数（ROTATELEFT/FF/GG/HH/II）和之前的FGHI一样，都是需要你进行SIMD并行化的
// 但是你需要注意的是#define的功能及其效果，可以发现这里的FGHI是没有返回值的，为什么呢？你可以查询#define的含义和用法

//--------串行部分-------//
#define ROTATELEFT(num, n) (((num) << (n)) | ((num) >> (32-(n))))


#define FF(a, b, c, d, x, s, ac) { \
  (a) += F ((b), (c), (d)) + (x) + ac; \
  (a) = ROTATELEFT ((a), (s)); \
  (a) += (b); \
}

#define GG(a, b, c, d, x, s, ac) { \
  (a) += G ((b), (c), (d)) + (x) + ac; \
  (a) = ROTATELEFT ((a), (s)); \
  (a) += (b); \
}
#define HH(a, b, c, d, x, s, ac) { \
  (a) += H ((b), (c), (d)) + (x) + ac; \
  (a) = ROTATELEFT ((a), (s)); \
  (a) += (b); \
}
#define II(a, b, c, d, x, s, ac) { \
  (a) += I ((b), (c), (d)) + (x) + ac; \
  (a) = ROTATELEFT ((a), (s)); \
  (a) += (b); \
}

//--------SIMD并行化部分-------//
/*
inline uint32x4_t ROTATELEFT_NEON(uint32x4_t vec, const int n) {
  uint32x4_t left = vshlq_n_u32(vec, n);          // 左移n位[1,6](@ref)
  uint32x4_t right = vshrq_n_u32(vec, 32 - n);   // 右移补位
  return vorrq_u32(left, right);                 // 按位或合并[1](@ref)
}

inline void FF_NEON(uint32x4_t &a, const uint32x4_t &b, 
  const uint32x4_t &c, const uint32x4_t &d,
  const uint32x4_t x, const int s, 
  const uint32_t ac) {
// 计算F(b,c,d) + x + ac
uint32x4_t F_val = F_neon(b, c, d);            // 调用NEON函数
uint32x4_t add_term = vaddq_u32(F_val, x);     // F + x
add_term = vaddq_u32(add_term, vdupq_n_u32(ac)); // 加常数
a = vaddq_u32(a, add_term);        // a += (...)
a = ROTATELEFT_NEON(a, s);         // 向量化循环左移
a = vaddq_u32(a, b);               // a += b
}

inline void GG_NEON(uint32x4_t &a, const uint32x4_t &b,
  const uint32x4_t &c, const uint32x4_t &d,
  const uint32x4_t x, const int s,
  const uint32_t ac) {
uint32x4_t G_val = G_neon(b, c, d);
uint32x4_t add_term = vaddq_u32(vaddq_u32(G_val, x), vdupq_n_u32(ac));
a = vaddq_u32(a, add_term);
a = ROTATELEFT_NEON(a, s);
a = vaddq_u32(a, b);
}

inline void HH_NEON(uint32x4_t &a,const uint32x4_t &b,
  const uint32x4_t &c, const uint32x4_t &d,
  const uint32x4_t x, const int s,
  const uint32_t ac) {
uint32x4_t H_val = H_neon(b, c, d);
uint32x4_t add_term = vaddq_u32(vaddq_u32(H_val, x), vdupq_n_u32(ac));
a = vaddq_u32(a, add_term);
a = ROTATELEFT_NEON(a, s);
a = vaddq_u32(a, b);
}

inline void II_NEON(uint32x4_t &a, const uint32x4_t &b,
  const uint32x4_t &c, const uint32x4_t &d,
  const uint32x4_t x, const int s,
  const uint32_t ac) {
uint32x4_t I_val = I_neon(b, c, d);
uint32x4_t add_term = vaddq_u32(vaddq_u32(I_val, x), vdupq_n_u32(ac));
a = vaddq_u32(a, add_term);
a = ROTATELEFT_NEON(a, s);
a = vaddq_u32(a, b);
}
*/

// 定义循环左移宏（兼容固定移位值）
#define ROTATELEFT_NEON(x, n) \
    vorrq_u32(vshlq_n_u32((x), (n)), vshrq_n_u32((x), 32 - (n)))

// 定义FF/GG/HH/II宏（兼容固定移位值）
#define FF_NEON(a, b, c, d, x, s, ac) { \
  uint32x4_t F_val = F_NEON((b), (c), (d)); \
  (a) = vaddq_u32((a), (F_val)); \
  (a) = vaddq_u32((a), (x)); \
  (a) = vaddq_u32((a), vdupq_n_u32(ac)); \
  (a) = ROTATELEFT_NEON((a), (s)); \
  (a) = vaddq_u32((a), (b)); \
}

#define GG_NEON(a, b, c, d, x, s, ac) { \
  uint32x4_t G_val = G_NEON((b), (c), (d)); \
  (a) = vaddq_u32((a), (G_val)); \
  (a) = vaddq_u32((a), (x)); \
  (a) = vaddq_u32((a), vdupq_n_u32(ac)); \
  (a) = ROTATELEFT_NEON((a), (s)); \
  (a) = vaddq_u32((a), (b)); \
}

#define HH_NEON(a, b, c, d, x, s, ac) { \
  uint32x4_t H_val = H_NEON((b), (c), (d)); \
  (a) = vaddq_u32((a), (H_val)); \
  (a) = vaddq_u32((a), (x)); \
  (a) = vaddq_u32((a), vdupq_n_u32(ac)); \
  (a) = ROTATELEFT_NEON((a), (s)); \
  (a) = vaddq_u32((a), (b)); \
}

#define II_NEON(a, b, c, d, x, s, ac) { \
  uint32x4_t I_val = I_NEON((b), (c), (d)); \
  (a) = vaddq_u32((a), (I_val)); \
  (a) = vaddq_u32((a), (x)); \
  (a) = vaddq_u32((a), vdupq_n_u32(ac)); \
  (a) = ROTATELEFT_NEON((a), (s)); \
  (a) = vaddq_u32((a), (b)); \
}

void MD5Hash(string input, bit32 *state);
void MD5HashFour(const string& pw0, const string& pw1, const string& pw2, const string& pw3, bit32 states[4][4]);