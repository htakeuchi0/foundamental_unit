/*
 * 実二次体K=Q(√m) (mは平方因子をもたない正整数) の基本単数を求めるプログラム．
 *
 * このモジュールは，
 * 実二次体K=Q(√m) (mは平方因子をもたない正整数) の基本単数を求めるプログラムである．
 *
 * [1]の§22.単数(2)より，K=Q(√m)について、判別式をDとすると、Kの単数はT, Uに関する
 * 以下の整数方程式
 *  
 *     T^2 - U^2D = ±4
 *
 * の解(T, U)を用いて，
 *
 *          T + U√D
 *     ε = ---------
 *             2
 *
 * と表される．特に最小正のU=uとそれに対応するT=tをとると，
 *
 *           t + u√D
 *     ε0 = ---------
 *              2
 *    
 * がKの基本単数(>1)となる．
 *
 * ただし，m≡2, 3 (mod 4) の場合は，
 *
 *     S^2 - U^2m = ±1    (*)
 *
 * の解(S, U)，最小解(s, u)を用いて，
 *
 *     ε = S + U√m,    ε0 = s + u√m
 *
 * と書ける．
 *
 * (*)の形の整数方程式はペル方程式と呼ばれ，
 * √mの連分数展開を用いた最小解の構成法が知られている．
 * 本モジュールでは，m≡2, 3 (mod 4) の場合はこの方法を用いて求める．
 *
 * m≡1 (mod 4) の場合でも，m = t^2 + 4 (tは正整数) と書ける場合，
 * 基本単数が簡単に計算できることが知られている [1, §22, 例3]. 
 * 本モジュールでは，この解法も採用して，実二次体の基本単数を求める．
 *
 * m≡1 (mod 4) でm = t^2 + 4 (tは正整数) と書けない場合でも
 * (1+√m)/2の連分数展開を用いて最小解が構成できることが知られている [2]．
 * 本稿では，[2]の結果を用いて，この場合の基本単数を求める．
 *
 * 参考文献:
 * [1] 石田 信，数学全書5 代数的整数論，森北出版株式会社，東京，1985.
 * [2] 有澤 健治，平方根の連分数とペル方程式 第3版，http://ar.nyx.link/cf/pell.pdf，2018
 *     (最終閲覧日: 2019/2/14).
 */

#ifndef FOUNDAMENTAL_UNIT_UNIT_H
#define FOUNDAMENTAL_UNIT_UNIT_H

#include <limits>
#include <iostream>
#include <stdexcept>

// GMPを使わない場合は次の行をコメントアウト
#define GMP

// 整数型の定義
#ifdef GMP
    #include <gmpxx.h>

    // 整数型
    using LongInteger = mpz_class;
    using SignedLongInteger = mpz_class;
#else
    // 整数型
    using LongInteger = unsigned long long int;
    using SignedLongInteger = long long int;
#endif  // #ifdef GMP

/** 配列の大きさ
 */
static const int ARRAY_SIZE = 100000;


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
bool SatisfiesRequirementsSquare(const LongInteger& a);


/** 与えられた整数が平方数ならtrueを返し、平方根も返す。
 *
 * @param a 平方数か判定したい整数
 * @param root 平方数の場合は平方根（そうでない場合は未定義）
 * @return 平方数ならtrue
 */
bool IsSquare(const LongInteger& a, LongInteger& root);


/** 実二次体K=Q(√m) (mは平方因子をもたない整数) の基本単数をナイーブな方法で計算する。
 *
 * @param m 実二次体K=Q(√m)のm
 * @param t 基本単数ε=(t+u√D)/2のt (ただし、DはKの判別式)
 * @param u 基本単数ε=(t+u√D)/2のu (ただし、DはKの判別式)
 * @return 基本単数のノルム
 */
int FoundamentalUnitNaive(int m, LongInteger& t, LongInteger& u);


/** 配列 a の a[start]〜a[end] が対称になっていればtrueを返す．
 *
 * @param a 配列
 * @param start 対称となる部分配列の初めの要素番号
 * @param start 対称となる部分配列の最後の要素番号
 * @return 配列が対称になっているときtrue
 */
bool IsCheckArray(int *a, int start, int end);


/** 平方根の整数部分を返す．
 * 
 * @param n 整数
 * @return nの平方根の整数部分
 */
int SquareRootIntegerPart(int n);


/** 整数mに対する(1+√m)/2の整数部分を返す．
 * 
 * @param m 整数
 * @return (1+√m)/2の整数部分
 */
int SquareRootIntegerPartExtended(int m);


/* 整数mに対する p_numer/denom + (q_numer/denom)*√mの整数部分を返す．
 * 
 * @param p_numer 整数
 * @param q_numer 整数
 * @param denom 整数
 * @param m 整数
 * @return p_numer/denom + (q_numer/denom)*√m の整数部分
 */
int SquareRootIntegerPartWithoutFloat(SignedLongInteger p_numer,
                                      SignedLongInteger q_numer,
                                      SignedLongInteger denom,
                                      int m);


/** 整数nに対する√nの連分数の係数を求め，循環節+1を返す．
 * 
 * @param n 整数
 * @param coeffs 連分数の係数
 * @param max_num_coeffs 連分数の係数の最大個数
 * @return 循環節+1の値（異常終了時は-1を返す）
 */
int ApproxContinuedFraction(int n, int *coeffs, int max_num_coeffs = ARRAY_SIZE);


/** 整数mに対する(1+√m)/2の連分数の係数を求め，循環節+1を返す．
 * 
 * @param m 整数
 * @param coeffs 連分数の係数
 * @param max_num_coeffs 連分数の係数の最大個数
 * @return 循環節+1の値（異常終了時は-1を返す）
 */
int ApproxContinuedFractionExtended(int m, int *coeffs, int max_num_coeffs = ARRAY_SIZE);


/** 連分数の係数から，分子と分母を計算して返す．
 *
 * @param coeffs 連分数の係数
 * @param len 係数の個数
 * @param numer 分子
 * @param denom 分母
 */
void CompContinuedFraction(int *coeffs, int len, 
                           SignedLongInteger& numer, SignedLongInteger& denom);


/** 実二次体K=Q(√m) (mは平方因子をもたない整数) の基本単数をペル方程式の解法を使って計算する。
 *
 * @param m 実二次体K=Q(√m)のm
 * @param t 基本単数ε=(t+u√D)/2のt (ただし、DはKの判別式)
 * @param u 基本単数ε=(t+u√D)/2のu (ただし、DはKの判別式)
 * @return 基本単数のノルム
 */
int FoundamentalUnitPellEq(int m, LongInteger& t, LongInteger& u);


/** 実二次体K=Q(√m) (mは平方因子をもたない正整数) の基本単数を
 * 拡張されたペル方程式の解法を使って計算する。
 *
 * @param m 実二次体K=Q(√m)のm
 * @param t 基本単数ε=(t+u√D)/2のt (ただし、DはKの判別式)
 * @param u 基本単数ε=(t+u√D)/2のu (ただし、DはKの判別式)
 * @return 基本単数のノルム
 */
int FoundamentalUnitPellEqExtended(int m, LongInteger& t, LongInteger& u);


/** 実二次体K=Q(√m)の基本単数を整えて表示する。
 *
 * @param m 実二次体K=Q(√m)のm
 * @param t 基本単数ε=(t+u√m)/2のt
 * @param u 基本単数ε=(t+u√m)/2のu
 * @param sign 基本単数のノルム
 * @param is_table_form 表形式で表示するときtrue
 */
void Show(int m, const LongInteger& t, const LongInteger& u, int sign, bool is_table_form);


/** 基本単数を表示する。
 *
 * @param m K=Q(√m) (mは平方因子をもたない) におけるm
 * @param is_table_form 表形式で表示するときtrue
 */
void DisplayFoundamentalUnit(int m);


/** 与えられた範囲の基本単数を表示する。
 *
 * @param min_num K=Q(√m) (mは平方因子をもたない) におけるmの最小値
 * @param max_num K=Q(√m) (mは平方因子をもたない) におけるmの最大値
 */
void DisplayFoundamentalUnits(int min_num = 2, int max_num = 200);

/** K=Q(√m)の基本単数(t+u√D)/2と表したときの，t, uを取得する
 *
 * @param m 実二次体K=Q(√m)のm
 * @param t 基本単数ε=(t+u√m)/2のt
 * @param u 基本単数ε=(t+u√m)/2のu
 * @return sign 基本単数のノルム
 */
int FoundamentalUnit(int m, LongInteger& t, LongInteger& u);

/** ヘルプを表示する
 */
void PrintUsage();

#endif // #ifndef FOUNDAMENTAL_UNIT_UNIT_H
