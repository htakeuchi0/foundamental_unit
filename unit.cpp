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
 *          T + U√D
 *     ε = ---------
 *             2
 *
 * と表される。特に最小正のU=uとそれに対応するT=tをとると、
 *
 *           t + u√D
 *     ε0 = ---------
 *              2
 *    
 * がKの基本単数(>1)となる。
 *
 * 参考文献:
 * [1] 石田 信，数学全書5 代数的整数論，森北出版株式会社，東京，1985.
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


/* 実二次体K=Q(√m) (mは平方因子をもたない整数) の判別式を返す。
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


/* 実二次体K=Q(√m) (mは平方因子をもたない整数) の基本単数を計算する。
 *
 * @param m 実二次体K=Q(√m)のm
 * @param t 基本単数ε=(t+u√D)/2のt (ただし、DはKの判別式)
 * @param u 基本単数ε=(t+u√D)/2のu (ただし、DはKの判別式)
 * @return 基本単数のノルム
 */
int FoundamentalUnit(int m, LongInteger& t, LongInteger& u) {
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


/* 実二次体K=Q(√m)の基本単数を整えて表示する。
 *
 * @param m 実二次体K=Q(√m)のm
 * @param t 基本単数ε=(t+u√D)/2のt (ただし、DはKの判別式)
 * @param u 基本単数ε=(t+u√D)/2のu (ただし、DはKの判別式)
 * @param d Kの判別式
 * @param sign 基本単数のノルム
 */
void Show(int m, const LongInteger& t, const LongInteger& u, 
        const LongInteger& d, int sign) {
    bool is_divisible_by_two = ((u & 0b1) == 0) && ((t & 0b1) == 0);

    // u, tがともに偶数の場合は約分した状態で表示する。
    if (is_divisible_by_two) {
        LongInteger t_divided = t >> 1;
        LongInteger u_divided = u >> 1;

        // uが1のときは1を省略する
        if (u_divided == 1) {
            std::cout << m << "," << sign << "," << t_divided << "+" << "√" << d << "\n";
        }
        else {
            std::cout << m << "," << sign << "," << t_divided << "+" 
                << u_divided << "√" << d << "\n";
        }
    }
    else {
        // uが1のときは1を省略する
        if (u == 1) {
            std::cout << m << "," << sign << ",(" << t << "+" << "√" << d << ")/2\n";
        }
        else {
            std::cout << m << "," << sign << ",(" << t << "+" << u << "√" << d << ")/2\n";
        }
    }
}


/* 基本単数を表示する。
 *
 * @param max_num K=Q(√m) (mは平方因子をもたない) におけるmの最大値
 */
void DisplayFoundamentalUnits(int max_num) {
    // ヘッダの表示
    std::cout << "# m, N(ε0), ε0\n";

    // 各K=Q(√m) (m=2,3,...,max_num) について基本単数を求める
    for (int m = 2; m <= max_num; m++) {
        // 平方部分をもつ場合は飛ばす
        if (SquarePart(m) != 1) {
            continue;
        }

        LongInteger t;
        LongInteger u;
        try {
            // 基本単数を求める
            int sign = FoundamentalUnit(m, t, u);

            // 判別式
            LongInteger d = Discriminant(m);

            // 判別式の平方部分
            LongInteger square_part_d = SquarePart(d);

            // u√dを整理
            d /= (square_part_d * square_part_d);
            u *= square_part_d;

            // 表示
            Show(m, t, u, d, sign);
        }
        catch (const std::exception& e) {
            // エラー内容を表示
            std::cout << m << ",0,(***" << e.what() << ")\n";
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
    int left = 1;
    int right = n;
    while (right - left > 1) {
        int middle = (left + right) >> 1;

        if (n < middle*middle) {
            right = middle;
        }
        else {
            left = middle;
        }
    }

    return left;
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
    // そこで，a[i]は整数であることを利用して，以下のようにω[i+1]を精度良く求める．
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
    //              α + β√n    αγ - βδn   βγ - αδ
    //    ω[i+1] = --------- = -------- + --------√n = p + q√n
    //              γ + δ√n    γ^2-δ^2n   γ^2-δ^2n
    //
    // を満たすp, qを浮動小数点数として求める．
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

    // a[1]以降を求める．
    for (int i = 1; i < max_num_coeffs; i++) {
        // 連分数計算
        long long int p0 = 1;
        long long int p1 = 0;
        long long int q0 = 0;
        long long int q1 = 1;

        for (int k = 1; k < i; k++) {
            long long int p2 = -coeffs[i - k] * p1 + p0;
            long long int q2 = -coeffs[i - k] * q1 + q0;

            p0 = p1;
            p1 = p2;
            q0 = q1;
            q1 = q2;
        }

        // α, β, γ, δの計算（すべて整数）
        long long int a = p0 - coeffs[0]*p1;
        long long int b = p1;
        long long int c = q0 - coeffs[0]*q1;
        long long int d = q1;

        // p, qの分母・分子を計算（整数）
        long long int p_numer = a*c - b*d*n;
        long long int q_numer = b*c - a*d;
        long long int denom = c*c - d*d*n;

        // p, qの計算（倍精度浮動小数点数）
        double p = static_cast<double>(p_numer) / denom;
        double q = static_cast<double>(q_numer) / denom;

        // Newton法．初期値はp + qnとする．
        double omega = p + q*n;
        int max_k = 100;
        for (int k = 0; k < max_k; k++) {
            omega -= (omega*omega - 2.0*p*omega + (p*p - q*q*n)) / (2.0*(omega - p));
        }

        // omegaの整数部分が連分数の係数になる
        coeffs[i] = static_cast<int>(omega);

        // 循環節が見つかったらループを終了
        if (coeffs[i] == 2*coeffs[0] && IsCheckArray(coeffs, 1, i-1)) {
            return i + 1;
        }
    }

    return -1;
}


/** 実験用メインメソッド
 */
int main() {
    // 基本単数を表示する
    DisplayFoundamentalUnits();

    return 0;
}
