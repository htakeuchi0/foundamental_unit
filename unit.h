/*
 * 二次体K=Q(√m) (mは平方因子をもたない) の基本単数を求めるプログラム。
 *
 * このモジュールは、[1]の§22.単数(2)に基づくナイーブな実装による、
 * 二次体K=Q(√m) (mは平方因子をもたない) の基本単数を求めるプログラムである。
 *
 * K=Q(√m)について、判別式をDとすると、Kの単数はT, Uに関する
 * 以下の整数方程式
 *  
 *     T^2 - U^2D = ±4
 *
 * の解(T, U)を用いて、
 *
 *          T - U√D
 *     ε = ---------
 *             2
 *
 * と表される。特に最小正のU=uとそれに対応するT=tをとると、
 *
 *           t - u√D
 *     ε0 = ---------
 *              2
 *    
 * がKの基本単数(>1)となる。
 *
 * 参考文献:
 * [1] 石田 信，数学全書5 代数的整数論，森北出版株式会社，東京，1985.
 */

#ifndef ALGEBRA_UNIT_H
#define ALGEBRA_UNIT_H

#include <iostream>
#include <stdexcept>

// GMPを使わない場合は次の行をコメントアウト
#define GMP

// 整数型の定義
#ifdef GMP
    #include <gmpxx.h>

    // 整数型
    using LongInteger = mpz_class;
#else
    // 整数型
    using LongInteger = unsigned long long int;
#endif  // #ifdef GMP


/** 静的配列のサイズを返す。
 *
 * @tparam TYPE 配列の型
 * @tparam SIZE 配列のサイズ
 * @param array TYPE型の静的配列
 * @return 配列の要素数
 */
template<typename TYPE, std::size_t SIZE>
std::size_t length(const TYPE(&)[SIZE]) {
    return SIZE;
}


/** 整数の平方部分を返す。
 *
 * @param a 整数
 * @param start 探索開始整数
 * @return a = bm^2 (bは平方数でない整数）における整数m
 */
LongInteger SquarePart(const LongInteger& a, const LongInteger& start = 2);



/** 実二次体K=Q(√m) (mは平方因子をもたない整数) の判別式を返す。
 *
 * @param m 実二次体K=Q(√m)におけるm
 * @return Kの判別式
 */
int Discriminant(int m);


/** 整数が平方数であるための必要条件を満たす場合trueを返す。
 *
 * 判定対象の整数の10進下2桁が、
 * 平方数の10進下2桁のパターンと合致する場合のみtrueを返す。
 *
 * trueが返された場合、それが平方数であるとは限らないが、
 * falseが返された場合、それは平方数ではないことがわかる。
 *
 * @param 平方数か判定したい整数
 * @param 必要条件を満たす場合true
 */
bool MeetRequirementsSquare(const LongInteger& a);


/** 与えられた整数が平方数ならtrueを返し、平方根も返す。
 *
 * @param a 平方数か判定したい整数
 * @param t 平方数の場合は平方根（そうでない場合は未定義）
 * @return 平方数ならtrue
 */
bool IsSquare(const LongInteger& a, LongInteger& t);


/** 実二次体K=Q(√m)の基本単数を整えて表示する。
 *
 * @param m 実二次体K=Q(√m)のm
 * @param u 基本単数ε=(t+u√D)/2のu (ただし、DはKの判別式)
 * @param t 基本単数ε=(t+u√D)/2のt (ただし、DはKの判別式)
 * @param d Kの判別式
 * @param sign 基本単数のノルム
 */
void Show(int m, const LongInteger& u, const LongInteger& t, const LongInteger& d, int sign);


/** 基本単数を表示する。
 *
 * @param max_num K=Q(√m) (mは平方因子をもたない) におけるmの最大値
 */
void DisplayFoundamentalUnits(int max_num = 200);

#endif // #ifndef ALGEBRA_UNIT_H
