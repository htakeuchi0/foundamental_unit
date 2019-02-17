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

#include "unit.h"


/* 整数の平方部分を返す。
 *
 * @param a 整数
 * @param start 探索開始整数
 * @return a = bm^2 (bは平方数でない整数）における整数m
 */
LongInteger SquarePart(const LongInteger& a, const LongInteger& start) {
    // 平方因子をもたない場合
    if (start * start > a) {
        return 1;
    }

    // 再帰的に探索
    for (LongInteger i = start; i < a; i++) {
        // 平方数で割り切れるか確認
        LongInteger j = i * i;
        if (a % j == 0) {
            return i * SquarePart(a / j, i);
        }
    }

    // 平方因子をもたない場合
    return 1;
}


/* 実二次体K=Q(√m) (mは平方因子をもたない正整数) の判別式を返す。
 *
 * @param m 実二次体K=Q(√m)におけるm
 * @return Kの判別式
 */
int Discriminant(int m) {
    return ((m & 0b11) == 1) ? m : (m << 2);
}


/* 整数が平方数であるための必要条件を満たす場合trueを返す。
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
bool SatisfiesRequirementsSquare(const LongInteger& a) {
    // 平方数の下2桁のパターン
    int lower_digits[] = {0, 1, 4, 9, 16, 21, 24, 25, 29, 36, 41, 44,
        49, 56, 61, 64, 69, 76, 81, 84, 89, 96};

    // 下2桁
    LongInteger a_lower_digit = a % 100;

    // 下2桁が平方数の下2桁リストにあるか二分法で検索する
    int left = 0;
    int right = length(lower_digits) - 1;
    LongInteger digit = -1;
    while (left <= right) {
        int middle = (left + right) >> 1;
        digit = lower_digits[middle];

        if (a_lower_digit == digit) {
            break;
        }
        else if (a_lower_digit < digit) {
            right = middle - 1;
        }
        else {
            left = middle + 1;
        }
    }

    // リストに含まれていれば、平方数の可能性がある。
    return (a_lower_digit == digit);
}


/* 与えられた整数が平方数ならtrueを返し、平方根も返す。
 *
 * @param a 平方数か判定したい整数
 * @param root 平方数の場合は平方根（そうでない場合は未定義）
 * @return 平方数ならtrue
 */
bool IsSquare(const LongInteger& a, LongInteger& root) {
    // 初期化
    root = -1;

    // 平方数の可能性がない場合はfalseを返して終了
    if (!SatisfiesRequirementsSquare(a)) {
        return false;
    }

    // 平方数候補の平方根を二分法で求める
    LongInteger left = 0;
    LongInteger right = a;
    while (left <= right) {
        LongInteger middle = (left + right) >> 1;

#ifndef GMP
        // オーバフロー判定
        if (middle >= (1LL << 32)) {
            throw std::runtime_error("Overflow has occurred while computing.");
        }
#endif

        LongInteger squared_middle = middle * middle;
        if (a == squared_middle) {
            root = middle;
            return true;
        }
        else if (a < squared_middle) {
            right = middle - 1;
        }
        else {
            left = middle + 1;
        }
    }

    return false;
}


/* 実二次体K=Q(√m) (mは平方因子をもたない正整数) の基本単数をナイーブな方法で計算する。
 *
 * @param m 実二次体K=Q(√m)のm
 * @param t 基本単数ε=(t+u√D)/2のt (ただし、DはKの判別式)
 * @param u 基本単数ε=(t+u√D)/2のu (ただし、DはKの判別式)
 * @return 基本単数のノルム
 */
int FoundamentalUnitNaive(int m, LongInteger& t, LongInteger& u) {
    // t, uの初期化
    t = -1;
    u = 1;

    // 判別式の計算
    int d = Discriminant(m);

    // 符号の初期化
    int sign = -1;

    // 方程式 T^2-U^2・D=±4を解く
    while (true) {
        // ±4+U^2・Dが負の場合
        if (sign == -1 && u*u*d < 4) {
            sign = -sign;
            continue;
        }

        // ±4+U^2・Dが平方数か判定し、平方数の場合、その平方根をtとする
        LongInteger squared_t = sign*4 + u*u*d;
        if (IsSquare(squared_t, t)) {
            return sign;
        }

        // ±4の符号の変換
        sign = -sign;

        // uを1増やす
        if (sign < 0) {
            u++;
        }
    }
}


/* 配列 a の a[start]〜a[end] が対称になっていればtrueを返す．
 *
 * @param a 配列
 * @param start 対称となる部分配列の初めの要素番号
 * @param start 対称となる部分配列の最後の要素番号
 */
bool IsCheckArray(int *a, int start, int end) {
    int i = start;
    int j = end;

    while (i < j) {
        if (a[i] != a[j]) {
            return false;
        }
        i++;
        j--;
    }

    return true;
}


/* 平方根の整数部分を返す．
 * 
 * @param n 整数
 * @return nの平方根の整数部分
 */
int SquareRootIntegerPart(int n) {
    LongInteger left = 1;
    LongInteger right = n;
    while (right - left > 1) {
        LongInteger middle = (left + right) >> 1;

        if (static_cast<unsigned int>(n) < middle*middle) {
            right = middle;
        }
        else {
            left = middle;
        }
    }

#ifdef GMP
    return left.get_si();
#else
    return left;
#endif // #ifdef GMP
}


/* 整数mに対する(1+√m)/2の整数部分を返す．
 * 
 * @param m 整数
 * @return (1+√m)/2の整数部分
 */
int SquareRootIntegerPartExtended(int m) {
    // m≧1より1≦(1+√m)/2≦m
    LongInteger left = 1;
    LongInteger right = m;
    while (right - left > 1) {
        LongInteger middle = (left + right) >> 1;

        LongInteger tmp = (middle << 1) - 1;
        if (static_cast<unsigned int>(m) < tmp*tmp) {
            right = middle;
        }
        else {
            left = middle;
        }
    }
#ifdef GMP
    return left.get_si();
#else
    return left;
#endif // #ifdef GMP
}

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
                                      int m) {
    LongInteger left = 1;
    LongInteger numer = p_numer + q_numer*m;
#ifdef GMP
    LongInteger right = static_cast<LongInteger>(numer.get_d()/denom);
#else
    LongInteger right = static_cast<LongInteger>(static_cast<double>(numer)/denom);
#endif // #ifdef GMP

    while (right - left > 1) {
        LongInteger middle = (left + right) >> 1;

        LongInteger lhs = denom*middle - p_numer;
        lhs *= lhs;
        LongInteger rhs = q_numer * q_numer * m;
        if (lhs > rhs) {
            right = middle;
        }
        else {
            left = middle;
        }
    }
#ifdef GMP
    return left.get_si();
#else
    return left;
#endif // #ifdef GMP
}


/* 整数nに対する√nの連分数の係数を求め，循環節+1を返す．
 * 
 * @param n 整数
 * @param coeffs 連分数の係数
 * @param max_num_coeffs 連分数の係数の最大個数
 * @return 循環節+1の値（異常終了時は-1を返す）
 */
int ApproxContinuedFraction(int n, int *coeffs, int max_num_coeffs) {
    // √nの連分数は，書き下すと，
    //                           1
    // √n = a[0] + -----------------------------
    //                              1 
    //             a[1] + ---------------------- 
    //                                1
    //                     a[2] + --------------
    //                               ...
    //                                 
    // となるため，以下の手続きで a[0], a[1], ... が求まる．
    //
    //   1: ω[0] := √n 
    //   2: i := 0
    //   3: while do
    //   4:     a[i] := floor(ω[i])
    //   5:     ω[i+1] := 1 / (ω[i] - a[i])
    //   6:     i := i + 1
    //   7: end while
    //
    // ただし，√nを浮動小数点数演算で求めたときに含まれる誤差や，
    // ω[i]->ω[i+1] の更新における計算誤差が，繰り返し計算で蓄積するため，
    // 上記手続きをそのまま実装すると正しく係数 a が求まらないことがある．
    // （作成者環境では√139で破綻した）
    //
    // そこで，a[i]は整数であることを利用して，以下のように
    // ω[i+1]の整数部分か，ω[i+1]そのものを精度良く求める．
    // (1) まず，以下の整数 α, β, γ, δ を求める．
    //
    //              α + β√n
    //    ω[i+1] = --------      (*)
    //              γ + δ√n
    //
    // ω[i+1]は連分数の形をとるため簡単に計算できる．
    // また，整数のみの演算で計算できるため，この計算中に誤差は含まれない．
    // 
    // (2) 次に，(*)の右辺を変形し，
    //    
    //              α + β√n    αγ - βδn   βγ - αδ      p_numer   q_numer
    //    ω[i+1] = --------- = -------- + --------√n = ------- + -------√n
    //              γ + δ√n    γ^2-δ^2n   γ^2-δ^2n      denom     denom
    //
    // を満たすp_numer, q_numer, denomを整数演算で求める．
    // 浮動小数点数演算を使わない場合は，このままω[i+1]の整数部分を求める．
    // したがって，浮動小数点数を使わない場合は，計算誤差は含まれない．
    //
    // 浮動小数点数演算を使う場合は, 
    // p = p_numer/denom, q = q_numer/denomを浮動小数点数として求める．
    //
    // (3) p + q√n をNewton法で求める．
    // 
    // Newton法は，近似値をωとすると，ω(new) := ω(old) + Δω の形となる．
    // 浮動小数点数演算では，Δω が ω(old)に比べて十分小さければ，
    // Δωの精度が悪くても，それが ω(new) の精度に大きく影響しないことが期待できる．
    //
    // p + q√n は2次方程式 ω^2 - 2pω + (p^2 - q^2n) = 0 の根なので，
    // qの符号に関わらず，初期値を p + qn としてNewton法を適用すればよい．
    // この方法では，平方根計算における誤差や，計算誤差の蓄積なく，
    // Newton法の更新式の性質とあわせ，精度良くω[i+1]が計算できることが期待される．
    // (ω[i+1]の整数部分が正しく求める程度の精度でよい点に注意されたい）
    
    // nの平方根の整数部分を計算
    int sqrt_int = SquareRootIntegerPart(n);
    coeffs[0] = static_cast<int>(sqrt_int);

    // 浮動小数点数演算を認めるか
#ifdef GMP
    bool use_float = false;
#else
    bool use_float = true;
#endif // #ifdef GMP


    // a[1]以降を求める．
    for (int i = 1; i < max_num_coeffs; i++) {
        // 連分数計算
        SignedLongInteger p0 = 1;
        SignedLongInteger p1 = 0;
        SignedLongInteger q0 = 0;
        SignedLongInteger q1 = 1;

        for (int k = 1; k < i; k++) {
            SignedLongInteger p2 = -coeffs[i - k] * p1 + p0;
            SignedLongInteger q2 = -coeffs[i - k] * q1 + q0;

            p0 = p1;
            p1 = p2;
            q0 = q1;
            q1 = q2;
        }

        // α, β, γ, δの計算（すべて整数）
        SignedLongInteger a = p0 - coeffs[0]*p1;
        SignedLongInteger b = p1;
        SignedLongInteger c = q0 - coeffs[0]*q1;
        SignedLongInteger d = q1;

        // p, qの分母・分子を計算（整数）
        SignedLongInteger p_numer = a*c - b*d*n;
        SignedLongInteger q_numer = b*c - a*d;
        SignedLongInteger denom = c*c - d*d*n;

        if (use_float) {
#ifdef GMP
            // p, qの計算（倍精度浮動小数点数）
            double p = p_numer.get_d() / denom.get_d();
            double q = q_numer.get_d() / denom.get_d();
#else
            // p, qの計算（倍精度浮動小数点数）
            double p = static_cast<double>(p_numer) / denom;
            double q = static_cast<double>(q_numer) / denom;
#endif // #ifdef GMP

            // Newton法．初期値はp + qnとする．
            double omega = p + q*n;
            int max_k = 100;
            for (int k = 0; k < max_k; k++) {
                omega -= (omega*omega - 2.0*p*omega + (p*p - q*q*n)) / (2.0*(omega - p));
            }

            // omegaの整数部分が連分数の係数になる
            coeffs[i] = static_cast<int>(omega);
        }
        else {
            coeffs[i] = SquareRootIntegerPartWithoutFloat(p_numer, q_numer, denom, n);
        }

        // 循環節が見つかったらループを終了
        if (coeffs[i] == 2*coeffs[0] && IsCheckArray(coeffs, 1, i-1)) {
            return i + 1;
        }
    }

    return -1;
}


/* 整数mに対する(1+√m)/2の連分数の係数を求め，循環節+1を返す．
 * 
 * @param m 整数
 * @param coeffs 連分数の係数
 * @param max_num_coeffs 連分数の係数の最大個数
 * @return 循環節+1の値（異常終了時は-1を返す）
 */
int ApproxContinuedFractionExtended(int m, int *coeffs, int max_num_coeffs) {
    // (1+√m)/2の連分数は，書き下すと，
    // 
    // 1 + √m                       1
    // ------ = a[0] + -----------------------------
    //   2                           1 
    //                 a[1] + ---------------------- 
    //                                1
    //                        a[2] + --------------
    //                                   ...
    //                                 
    // となるため，以下の手続きで a[0], a[1], ... が求まる．
    //
    //   1: ω[0] := (1+√m)/2 
    //   2: i := 0
    //   3: while do
    //   4:     a[i] := floor(ω[i])
    //   5:     ω[i+1] := 1 / (ω[i] - a[i])
    //   6:     i := i + 1
    //   7: end while
    //
    // ただし，√mを浮動小数点数演算で求めたときに含まれる誤差や，
    // ω[i]->ω[i+1] の更新における計算誤差が，繰り返し計算で蓄積するため，
    // 上記手続きをそのまま実装すると正しく係数 a が求まらないことがある．
    //
    // そこで，a[i]は整数であることを利用して，以下のように
    // ω[i+1]の整数部分か，ω[i+1]そのものを精度良く求める．
    // (1) まず，以下の整数 α, β, γ, δ を求める．
    //
    //              α + β√n
    //    ω[i+1] = --------      (*)
    //              γ + δ√n
    //
    // ω[i+1]は連分数の形をとるため簡単に計算できる．
    // また，整数のみの演算で計算できるため，この計算中に誤差は含まれない．
    // 
    // (2) 次に，(*)の右辺を変形し，
    //    
    //              α + β√n    αγ - βδn   βγ - αδ      p_numer   q_numer
    //    ω[i+1] = --------- = -------- + --------√n = ------- + -------√n
    //              γ + δ√n    γ^2-δ^2n   γ^2-δ^2n      denom     denom
    //
    // を満たすp_numer, q_numer, denomを整数演算で求める．
    // 浮動小数点数演算を使わない場合は，このままω[i+1]の整数部分を求める．
    // したがって，浮動小数点数を使わない場合は，計算誤差は含まれない．
    //
    // 浮動小数点数演算を使う場合は, 
    // p = p_numer/denom, q = q_numer/denomを浮動小数点数として求める．
    //
    // (3) p + q√n をNewton法で求める．
    // 
    // Newton法は，近似値をωとすると，ω(new) := ω(old) + Δω の形となる．
    // 浮動小数点数演算では，Δω が ω(old)に比べて十分小さければ，
    // Δωの精度が悪くても，それが ω(new) の精度に大きく影響しないことが期待できる．
    //
    // p + q√n は2次方程式 ω^2 - 2pω + (p^2 - q^2n) = 0 の根なので，
    // qの符号に関わらず，初期値を p + qn としてNewton法を適用すればよい．
    // この方法では，平方根計算における誤差や，計算誤差の蓄積なく，
    // Newton法の更新式の性質とあわせ，精度良くω[i+1]が計算できることが期待される．
    // (ω[i+1]の整数部分が正しく求める程度の精度でよい点に注意されたい）
    
    // nの平方根の整数部分を計算
    int sqrt_int = SquareRootIntegerPartExtended(m);
    coeffs[0] = static_cast<int>(sqrt_int);

    // 浮動小数点数演算を認めるか
#ifdef GMP
    bool use_float = false;
#else
    bool use_float = true;
#endif // #ifdef GMP

    // a[1]以降を求める．
    for (int i = 1; i < max_num_coeffs; i++) {
        // 連分数計算
        SignedLongInteger p0 = 1;
        SignedLongInteger p1 = 0;
        SignedLongInteger q0 = 0;
        SignedLongInteger q1 = 1;

        for (int k = 1; k < i; k++) {
            SignedLongInteger p2 = -coeffs[i - k] * p1 + p0;
            SignedLongInteger q2 = -coeffs[i - k] * q1 + q0;

            p0 = p1;
            p1 = p2;
            q0 = q1;
            q1 = q2;
        }

        // α, β, γ, δの計算（すべて整数）
        SignedLongInteger tmp = (coeffs[0] << 1) - 1;
        SignedLongInteger a = (p0 << 1) - tmp*p1;
        SignedLongInteger b = p1;
        SignedLongInteger c = (q0 << 1) - tmp*q1;
        SignedLongInteger d = q1;

        // p, qの分母・分子を計算（整数）
        SignedLongInteger p_numer = a*c - b*d*m;
        SignedLongInteger q_numer = b*c - a*d;
        SignedLongInteger denom = c*c - d*d*m;

        if (use_float) {
#ifdef GMP
            // p, qの計算（倍精度浮動小数点数）
            double p = p_numer.get_d() / denom.get_d();
            double q = q_numer.get_d() / denom.get_d();
#else
            // p, qの計算（倍精度浮動小数点数）
            double p = static_cast<double>(p_numer) / denom;
            double q = static_cast<double>(q_numer) / denom;
#endif // #ifdef GMP

            // Newton法．初期値はp + qmとする．
            double omega = p + q*m;
            int max_k = 100;
            for (int k = 0; k < max_k; k++) {
                omega -= (omega*omega - 2.0*p*omega + (p*p - q*q*m)) / (2.0*(omega - p));
            }

            // omegaの整数部分が連分数の係数になる
            coeffs[i] = static_cast<int>(omega);
        }
        else {
            coeffs[i] = SquareRootIntegerPartWithoutFloat(p_numer, q_numer, denom, m);
        }

        // 循環節が見つかったらループを終了
        if (coeffs[i] == 2*coeffs[0] - 1 && IsCheckArray(coeffs, 1, i-1)) {
            return i + 1;
        }
    }

    return -1;
}


/* 連分数の係数から，分子と分母を計算して返す．
 *
 * @param coeffs 連分数の係数
 * @param len 係数の個数
 * @param numer 分子
 * @param denom 分母
 */
void CompContinuedFraction(int *coeffs, int len, 
                           SignedLongInteger& numer, SignedLongInteger& denom) {
    // 連分数計算
    // ユークリッドの互除法の逆算で計算する．
    SignedLongInteger p0 = 1;
    SignedLongInteger p1 = coeffs[0];
    SignedLongInteger q0 = 0;
    SignedLongInteger q1 = 1;

    for (int i = 1; i < len; i++) {
        SignedLongInteger p2 = coeffs[i] * p1 + p0;
        SignedLongInteger q2 = coeffs[i] * q1 + q0;

        p0 = p1;
        p1 = p2;
        q0 = q1;
        q1 = q2;
    }

    numer = p1;
    denom = q1;
}


/* 実二次体K=Q(√m) (mは平方因子をもたない正整数) の基本単数をペル方程式の解法を使って計算する。
 *
 * @param m 実二次体K=Q(√m)のm
 * @param t 基本単数ε=(t+u√D)/2のt (ただし、DはKの判別式)
 * @param u 基本単数ε=(t+u√D)/2のu (ただし、DはKの判別式)
 * @return 基本単数のノルム
 */
int FoundamentalUnitPellEq(int m, LongInteger& t, LongInteger& u) {
    // m≡2, 3 (mod 4) のとき，K=Q(√m)の基本単数ε0は，
    // 整数方程式 S^2 - U^2m = ±1 の最小整数解 (s, u) を計算し，
    // ε0 = s + u√m として求められる．
    //
    // 方程式 S^2 - U^2m = ±1 はペル方程式と呼ばれ，以下のように
    // √m の連分数展開を用いて最小解を計算できることが知られている．
    //
    // ペル方程式の最小解の構成法：
    //   √m = [a0; (a1, a2, ..., al)] (括弧内は循環節を表す）とするとき，
    //   p/q = [a0; a1, a2, ..., a(l-1)] を既約分数とすると，
    //   (s, u) = (p, q)が解となる．
    //
    // 連分数：
    //   [a0; a1, ..., al, ...]
    //
    //                        1
    //      := a0 + ---------------------
    //                          1 
    //               a1 + ---------------
    //                     ...
    //                                1
    //                         al + -----
    //                               ...
    SignedLongInteger p;
    SignedLongInteger q;
    
    // √mの連分数展開を計算する
    int coeffs[ARRAY_SIZE];
    int len = ApproxContinuedFraction(m, coeffs);
    if (len == -1) {
        throw std::runtime_error("*** [ERROR] Could not compute. (Computation might be possible by setting ARRAY_SIZE to a larger value.)");
    }

    // [a0; a1, ..., a(l-1)]を計算する
    CompContinuedFraction(coeffs, len-1, p, q); 

    // 表示処理統一のために T^2 - U^2D = ±4 の形式にしておく
    t = p << 1;
    u = q << 1;
    LongInteger sign_long_int = p*p - q*q*m;

#ifdef GMP
    int sign = sign_long_int.get_si();
#else
    int sign = sign_long_int;
#endif // #ifdef GMP

    return sign;
}


/* 実二次体K=Q(√m) (mは平方因子をもたない正整数) の基本単数を
 * 拡張されたペル方程式の解法を使って計算する。
 *
 * @param m 実二次体K=Q(√m)のm
 * @param t 基本単数ε=(t+u√D)/2のt (ただし、DはKの判別式)
 * @param u 基本単数ε=(t+u√D)/2のu (ただし、DはKの判別式)
 * @return 基本単数のノルム
 */
int FoundamentalUnitPellEqExtended(int m, LongInteger& t, LongInteger& u) {
    // m≡1 (mod 4) のとき，K=Q(√m)の基本単数ε0は，
    // 整数方程式 T^2 - U^2m = ±4 の最小整数解 (t, u) を計算し，
    // ε0 = (t + u√m)/2 として求められる．
    //
    // 方程式 T^2 - U^2m = ±4 も右辺が±1と同様ペル方程式と呼ばれるが，
    // 通常の右辺が±1の場合と全く同じ方法ではこの方程式は解けない．
    // ただし，(1+√m)/2 の連分数展開を用いて最小解を計算できることが，
    // 有澤によって主張されている．
    // (http://ar.nyx.link/cf/pell.pdf)
    //
    // 右辺が±4となるペル方程式（拡張ペル方程式）の最小解の構成法：
    //   (1+√m)/2 = [a0; (a1, a2, ..., al)] (括弧内は循環節を表す）とするとき，
    //   p/q = [a0; a1, a2, ..., a(l-1)] を既約分数とすると，
    //   (t, u) = (2p-q, q)が解となる．
    //
    // 連分数：
    //   [a0; a1, ..., al, ...]
    //
    //                        1
    //      := a0 + ---------------------
    //                          1 
    //               a1 + ---------------
    //                     ...
    //                                1
    //                         al + -----
    //                               ...
    SignedLongInteger p;
    SignedLongInteger q;
    
    // (1+√m)/2の連分数展開を計算する
    int coeffs[ARRAY_SIZE];
    int len = ApproxContinuedFractionExtended(m, coeffs);
    if (len == -1) {
        throw std::runtime_error("*** [ERROR] Could not compute. (Computation might be possible by setting ARRAY_SIZE to a larger value.)");
    }

    // [a0; a1, ..., a(l-1)]を計算する
    CompContinuedFraction(coeffs, len-1, p, q); 

    // T^2 - U^2D = ±4 のT, Uを求める
    t = (p << 1) - q;
    u = q;
    LongInteger sign_long_int = (t*t - u*u*m) >> 2;

#ifdef GMP
    int sign = sign_long_int.get_si();
#else
    int sign = sign_long_int;
#endif // #ifdef GMP

    return sign;
}


/* 実二次体K=Q(√m)の基本単数を整えて表示する。
 *
 * @param m 実二次体K=Q(√m)のm
 * @param t 基本単数ε=(t+u√m)/2のt
 * @param u 基本単数ε=(t+u√m)/2のu
 * @param sign 基本単数のノルム
 * @param is_table_form 表形式で表示するときtrue
 */
void Show(int m, const LongInteger& t, const LongInteger& u, int sign, bool is_table_form) {
    bool is_divisible_by_two = ((u & 0b1) == 0) && ((t & 0b1) == 0);

    // u, tがともに偶数の場合は約分した状態で表示する。
    if (is_divisible_by_two) {
        LongInteger t_divided = t >> 1;
        LongInteger u_divided = u >> 1;

        // uが1のときは1を省略する
        if (u_divided == 1) {
            if (is_table_form) {
                std::cout << m << "," << sign << "," << t_divided << "+" << "√" << m << "\n";
            }
            else {
                std::cout << "K:       " << "Q(√" << m << ")\n";
                std::cout << "ε0:      " << t_divided << "+" << "√" << m << "\n";
                std::cout << "N_K(ε0): " << sign << "\n";
            }
        }
        else {
            if (is_table_form) {
                std::cout << m << "," << sign << "," << t_divided << "+" 
                    << u_divided << "√" << m << "\n";
            }
            else {
                std::cout << "K:       " << "Q(√" << m << ")\n";
                std::cout << "ε0:      " << t_divided << "+" << u_divided << "√" << m << "\n";
                std::cout << "N_K(ε0): " << sign << "\n";
            }
        }
    }
    else {
        // uが1のときは1を省略する
        if (u == 1) {
            if (is_table_form) {
                std::cout << m << "," << sign << ",(" << t << "+" << "√" << m << ")/2\n";
            }
            else {
                std::cout << "K:       " << "Q(√" << m << ")\n";
                std::cout << "ε0:      " << "(" << t << "+" << "√" << m << ")/2\n";
                std::cout << "N_K(ε0): " << sign << "\n";
            }
        }
        else {
            if (is_table_form) {
                std::cout << m << "," << sign << ",(" << t << "+" << u << "√" << m << ")/2\n";
            }
            else {
                std::cout << "K:       " << "Q(√" << m << ")\n";
                std::cout << "ε0:      " << "(" << t << "+" << u << "√" << m << ")/2\n";
                std::cout << "N_K(ε0): " << sign << "\n";
            }
        }
    }
}


/** 基本単数を表示する。
 *
 * @param m K=Q(√m) (mは平方因子をもたない) におけるm
 * @param is_table_form 表形式で表示するときtrue
 */
void DisplayFoundamentalUnit(int m, bool is_table_form) {
    LongInteger t;
    LongInteger u;
    try {
        // 基本単数を求める
        int sign = FoundamentalUnit(m, t, u);

        // 表示
        Show(m, t, u, sign, is_table_form);
    }
    catch (const std::exception& e) {
        // エラー内容を表示
        if (is_table_form) {
            std::cout << m << ",0,(***" << e.what() << ")\n";
        }
        else {
            std::cout << "K:     " << "Q(√" << m << ")\n";
            std::cout << "*** [ERROR] " << e.what() << "\n";
        }
    }
}


/** 与えられた範囲の基本単数を表示する。
 *
 * @param min_num K=Q(√m) (mは平方因子をもたない) におけるmの最小値
 * @param max_num K=Q(√m) (mは平方因子をもたない) におけるmの最大値
 */
void DisplayFoundamentalUnits(int min_num, int max_num) {
    // ヘッダの表示
    if (min_num <= max_num) {
        std::cout << "# m, N_K(ε0), ε0\n";
    }

    // 各K=Q(√m) (m=2,3,...,max_num) について基本単数を求める
    for (int m = min_num; m <= max_num; m++) {
        // 平方部分をもつ場合は飛ばす
        if (SquarePart(m) != 1 || m == 1) {
            continue;
        }

        DisplayFoundamentalUnit(m, true);
    }
}


/* K=Q(√m)の基本単数(t+u√D)/2と表したときの，t, uを取得する
 *
 * @param m 実二次体K=Q(√m)のm
 * @param t 基本単数ε=(t+u√m)/2のt
 * @param u 基本単数ε=(t+u√m)/2のu
 * @return sign 基本単数のノルム
 */
int FoundamentalUnit(int m, LongInteger& t, LongInteger& u) {
#ifdef GMP
    int max_m = std::numeric_limits<int>::max();
#else
    int max_m = 523;
#endif  // #ifdef GMP

    // 負数判定
    if (m <= 0) {
        throw std::runtime_error("The value of m must be positive.");
    }

    // オーバフロー判定
    if (m > max_m) {
        throw std::runtime_error("The value of m must be " + std::to_string(max_m) + " or less.");
    }

    // 平方部分をもつか判定
    if (SquarePart(m) != 1) {
        throw std::runtime_error("The value of m has a square part.");
    }

    int sign = 0;

    // 基本単数を求める
    if ((m & 0b10) != 0) {
        // m≡2, 3 (mod 4) のとき，K=Q(√m)の基本単数ε0は，
        // 整数方程式 S^2 - U^2m = ±1 の最小整数解 (s, u) を計算し，
        // ε0 = s + u√m として求められる．
        //
        // 方程式 S^2 - U^2m = ±1 はペル方程式と呼ばれ，以下のように
        // √m の連分数展開を用いて最小解を計算できることが知られている．
        sign = FoundamentalUnitPellEq(m, t, u);
    }
    else {
        // m≡1 (mod 4) のとき，K=Q(√m)の基本単数ε0は，
        // 整数方程式 T^2 - U^2D = ±4 の最小整数解 (t, u) を計算し，
        // 基本単数 ε0 = (t + u√D)/2 を求められる．
        // ただし，m = t^2 + 4と書ける場合は，より簡単に計算できる．
        if (IsSquare(m - 4, t)) {
            // m = t^2 + 4の場合，ε0 = (t + √m)/2 となり，
            // N_K(ε0) = (t^2 - m)/4 = -1 である．
            u = 1;
            sign = -1;
        }
        else {
            // m≡1 (mod 4) のとき，K=Q(√m)の基本単数ε0は，
            // 整数方程式 T^2 - U^2m = ±4 の最小整数解 (t, u) を計算し，
            // ε0 = (t + u√m)/2 として求められる．
            //
            // 方程式 T^2 - U^2m = ±4 はペル方程式と呼ばれ，以下のように
            // (1+√m)/2 の連分数展開を用いて最小解を計算できることが知られている．
            sign = FoundamentalUnitPellEqExtended(m, t, u);
        }
    }

    return sign;
}

/* ヘルプを表示する
 */
void PrintUsage() {
    std::cout << "Usage:\n";
    std::cout << "  ./unit -h               show this help message\n";
    std::cout << "  ./unit --help           show this help message\n";
    std::cout << "\n";
    std::cout << "  ./unit                  display foundamental units of Q(√m)"
                 " for m = [2..199]\n";
    std::cout << "\n";
    std::cout << "  ./unit m                display a foundamental unit of Q(√m)\n";
    std::cout << "\n";
    std::cout << "  ./unit min_m max        display foundamental units of Q(√m)"
                 " for m = [min_m..max_m]\n";
    std::cout << "\n";
    std::cout << "Examples:\n";
    std::cout << "  ./unit 2\n";
    std::cout << "  ./unit 2 200\n";
}

/** 実験用メインメソッド
 */
int main(int argc, char *argv[]) {
    int m = 0;
    int min_m = 0;
    int max_m = 0;
    int int_max = std::numeric_limits<int>::max();

    // ヘルプの表示
    if (argc > 1) {
        std::string arg_str = argv[1];
        if (argv[1] == std::string("--help") || argv[1] == std::string("-h")) {
            PrintUsage();
            return 0;
        }
    }

    switch (argc) {
    case 1:
        // 基本単数のリストを表示する
        DisplayFoundamentalUnits();
        break;

    case 2:
        // intに丸める前にオーバーフロー判定
        if (std::atoll(argv[1]) > int_max) {
            std::cout << "*** [ERROR] At least the value of m must be INT_MAX or less.\n";
            return -1;
        }

        // 基本単数を表示する
        m = std::atoi(argv[1]);
        DisplayFoundamentalUnit(m, false);
        break;

    default:
        // intに丸める前にオーバーフロー判定
        if (std::atoll(argv[1]) > int_max || std::atoll(argv[2]) > int_max) {
            std::cout << "*** [ERROR] At least the value of m must be INT_MAX or less.\n";
            return -1;
        }

        // 基本単数を表示する
        min_m = std::atoi(argv[1]);
        max_m = std::atoi(argv[2]);
        DisplayFoundamentalUnits(min_m, max_m);
        break;
    }

    return 0;
}
