/**
Boost::uBlas目指そうと思ったけどuBlasそんなに使ったことないし、
じゃあ独自的な方向でいこうって感じのExpression Templateを使った行列ライブラリ。

特徴は
・ET使ってるから遅延評価してるらしい
・ET使ってるからコンパイラ酷使するらしい
・静的大きさがメインだけど、もちろん動的大きさも視野に入れて設計したらしい
・identityとかonesとか、特殊な行列はそもそも大きさを持たないらしい
・LU分解とそれを使った逆行列ぐらいはそのうち入れるらしい
・疎行列とかそういう特殊行列も入れたいらしい
・そんな時間ない気がするらしい
・とにかく、現時点ではコア機能しか使えないらしい
・開発動機は「面白そう」。
・実行速度ってどうしたら速くなるのか開発者はしらないらしい
    開発者は趣味でプログラミングしてるから、本気で数値計算を勉強したこと無いらしい
・素直にBLAS, LAPACKをつかいましょう。

TODO:
外積,
LU分解, 逆行列,
共役転置

型としての行列の種別
    帯行列, 対称行列, エルミート行列, 疎行列,

高速化とか
・Blas

特定行列の最適化(identity * A == AとかA.transpose.transpose == A)


Author: Kazuki Komatsu

License: NYSL
*/

module gigue.core;

import std.algorithm,
       std.array,
       std.conv,
       std.format,
       std.functional,
       std.range,
       std.traits,
       std.string;

version(unittest) import std.stdio;


/**
デフォルトでの行列が列優先か行優先のどちらでメモリ上に配置されるかをDlinearColumnMajorという
*/
/*
version(ExtmlColumnMajor)
    enum Major defaultMajor = Major.column;
else
    enum Major defaultMajor = Major.row;
*/

//version(unittest) void main(){}

enum size_t wild = 0;       /// inferableMatrix

/**
四則演算が定義されている型。
inout(ubyte)からinout(creal)などまで。
*/
template isScalar(T){
    enum bool isScalar = is(typeof(
        {
            T* a;

            {auto _t1 = *a + *a;}
            {auto _t1 = *a - *a;}
            {auto _t1 = *a * *a;}
            {auto _t1 = *a / *a;}

            //bool b = a < 0;
            bool b = *a == *a;
        }));
}
unittest{
    import std.bigint, std.numeric, std.typetuple;

    alias TT = TypeTuple!(ubyte, ushort, uint, ulong,
                           byte,  short,  int,  long,
                          float, double, real,
                          creal, cfloat, cdouble/*,*/
                          /*BigInt*/);

    foreach(T; TT)
    {
        static assert(isScalar!(T));
        static assert(isScalar!(const(T)));
        static assert(isScalar!(immutable(T)));
    }

    static assert(isScalar!BigInt);

    static assert(!isScalar!(CustomFloat!16));
}


/**
行列の格納方式
*/
enum Major{
    row,
    column,
}

enum defaultMajor = Major.row;


/**
行列型か判定します。
行列型Tは次のコードがコンパイル可能です。
*/
template isMatrix(T)
{
    enum bool isMatrix = is(typeof({
            T m;

            size_t rsize = m.rows;
            size_t csize = m.cols;

            auto e = m[0, 0];
            //static assert(isScalar!(typeof(e)));
        }));
}
unittest{
    static struct S
    {
        enum rows = 1;
        enum cols = 1;

        auto opIndex(size_t i, size_t j)
        {
            return i + j;
        }
    }

    static assert(isMatrix!S);
}
unittest{
    static struct S
    {
        enum rows = 0;
        enum cols = 0;

        auto opIndex(size_t i, size_t j)
        {
            return i + j;
        }
    }

    static assert(isMatrix!S);
}
unittest{
    import std.bigint, std.typetuple;
    alias TT = TypeTuple!(ubyte, ushort, uint, ulong,
                           byte,  short,  int,  long,
                          float, double, real,
                          creal, cfloat, cdouble/*,*/
                          /*BigInt*/);

    static struct M(T, size_t r, size_t c)
    {
        enum rows = r;
        enum cols = c;

        auto opIndex(size_t i, size_t j){return T.init;}
    }

    foreach(T; TT)
    {
        static assert(isMatrix!(M!(T, 3, 3)));
        static assert(isMatrix!(M!(const(T), 3, 3)));
        static assert(isMatrix!(M!(immutable(T), 3, 3)));
    }
}


/**
行のサイズが静的に決定されるか
*/
template hasStaticRows(T)
{
    enum bool hasStaticRows = is(typeof({
            enum rsize = T.rows;
            static assert(rsize > 0);
        }));
}


/**
列のサイズが静的に決定されるか
*/
template hasStaticColumns(T)
{
    enum bool hasStaticColumns = is(typeof({
            enum csize = T.cols;
            static assert(csize > 0);
        }));
}


/**
サイズが静的に決定されるか
*/
template hasStaticLength(T)
{
    enum bool hasStaticLength = is(typeof({
            enum size = T.length;
            static assert(size > 0);
        }));
}


/**
行のサイズが動的に変化する
*/
template hasDynamicRows(T)
{
    enum bool hasDynamicRows = !is(typeof({
            enum rsize = T.rows;
        }));
}


/**
列のサイズが動的に変化する
*/
template hasDynamicColumns(T)
{
    enum bool hasDynamicColumns = !is(typeof({
            enum csize = T.cols;
        }));
}


/**
サイズが動的に変化する
*/
template hasDynamicLength(T)
{
    enum bool hasDynamicLength = !is(typeof({
            enum size = T.length;
        }));
}


/**

*/
struct InferredResult
{
    bool isValid;
    size_t rows;
    size_t cols;
}


/**
サイズが推論可能
*/
template isInferableMatrix(T)
{
    enum bool isInferableMatrix = isMatrix!T && is(typeof({
            static assert(T.rows == wild);
            static assert(T.cols == wild);

            enum InferredResult result = T.inferSize(wild, wild);

            static if(result.isValid)
            {
                enum inferredRows = result.rows;
                enum inferredCols = result.cols;
            }

            enum someValue = 4; //非負の整数
            static assert(!T.inferSize(wild, wild).isValid);
            static assert(T.inferSize(wild, someValue).isValid);
            static assert(T.inferSize(wild, someValue).isValid);
        }));
}


/**
ベクトル型かどうか判定します。
*/
template isVector(V)
{
    enum isVector = isMatrix!V && !isInferableMatrix!V
                 && is(typeof({
                        static if(hasStaticRows!V)
                        {
                            static if(hasStaticColumns!V)
                                static assert(V.rows == 1 || V.cols == 1);
                            else
                                static assert(V.rows == 1);
                        }
                        else
                            static assert(V.cols == 1);

                        V v;
                        size_t size = v.length;

                        auto a = v[0, 0];
                        auto b = v[0];
                        static assert(is(typeof(a) == typeof(b)));
                    }));
}
unittest{
    static struct V
    {
        enum rows = 1;
        enum cols = 3;
        enum length = 3;

        auto opIndex(size_t i, size_t j){return j;}
        auto opIndex(size_t i){return i;}
    }

    static assert(isVector!V);
}


/**
行列型の要素の型を取得します。

Example:
---
static assert(is(ElementType!(NMatrix2S!int) == int));
---
*/
template ElementType(A) if(isMatrix!A)
{
    alias typeof(A.init[0, 0]) ElementType;
}
unittest{
    static struct S
    {
        enum rows = 1;
        enum cols = 1;

        int opIndex(size_t, size_t)
        {
            return 1;
        }
    }

    static assert(is(ElementType!S == int));
}


/**
その行列の要素がstd.algorithm.swapを呼べるかどうかをチェックします。

Example:
---
static assert(hasLvalueElements!(Matrix3i));
---
*/
template hasLvalueElements(A)if(isMatrix!A)
{
    enum hasLvalueElements = is(typeof({
            import std.algorithm : swap;

            A a;
            swap(a[0, 0], a[0, 0]);
        }));
}
unittest{
    static struct M
    {
        enum rows = 1;
        enum cols = 1;

        ref int opIndex(size_t i, size_t j)
        {
            static int a;
            return a;
        }
    }

    static assert(hasLvalueElements!(M));
}


/**
A型の行列の要素にElementType!A型の値が代入可能かどうかをチェックします。

Example:
---
static assert(hasAssignableElements!(Matrix3i));
---
*/
template hasAssignableElements(A)if(isMatrix!A)
{
    enum hasAssignableElements = is(typeof({
            ElementType!A e;
            A a;
            a[0, 0] = e;
        }));
}
unittest{
    static struct M
    {
        enum rows = 1;
        enum cols = 1;

        auto opIndex(size_t i, size_t j)
        {
            return i;
        }

        void opIndexAssign(typeof(this[0, 0]) a, size_t i, size_t j){}
    }

    static assert(hasAssignableElements!(M));
}


/**
2つの行列が演算可能かどうかを判定するのに使います。
aとbが等しいか、もしくはどちらかが0であるとtrueとなります。

Example:
---
static assert(isEqOrEitherEq0(0, 0));
static assert(isEqOrEitherEq0(1, 1));
static assert(isEqOrEitherEq0(0, 1));
static assert(isEqOrEitherEq0(1, 0));
---
*/
private bool isEqOrEitherEq0(alias pred = "a == b", T)(T a, T b)
{
    return binaryFun!pred(a, b) || a == 0 || b == 0;
}
unittest{
    static assert(isEqOrEitherEq0(0, 0));
    static assert(isEqOrEitherEq0(1, 1));
    static assert(isEqOrEitherEq0(0, 1));
    static assert(isEqOrEitherEq0(1, 0));
}



/**

*/
template isValidOperator(L, string op, R)
{
    static if(isMatrix!L && isMatrix!R)
        enum isValidOperator = is(typeof(mixin("L.init[0, 0] " ~ op ~ " R.init[0, 0]"))) && isValidOperatorImpl!(L, op, R);
    else static if(isMatrix!L)
        enum isValidOperator = is(typeof(mixin("L.init[0, 0] " ~ op ~ " R.init"))) && isValidOperatorImpl!(L, op, R);
    else static if(isMatrix!R)
        enum isValidOperator = is(typeof(mixin("L.init " ~ op ~ " R.init[0, 0]"))) && isValidOperatorImpl!(L, op, R);
}


template isValidOperatorImpl(L, string op, R)
if(isMatrix!L && isMatrix!R && op != "*")
{
    static if(op != "+" && op != "-")
        enum isValidOperatorImpl = false;
    else static if(isInferableMatrix!L && isInferableMatrix!R)
        enum isValidOperatorImpl = true;
    else static if(isInferableMatrix!L)
    {
        static if(hasStaticRows!R)
            enum isValidOperatorImpl = isValidOperatorImpl!(Inferred!(L, R.rows, wild), op, R);
        else static if(hasStaticColumns!R)
            enum isValidOperatorImpl = isValidOperatorImpl!(Inferred!(L, wild, R.cols), op, R);
        else
            enum isValidOperatorImpl = true;
    }
    else static if(isInferableMatrix!R)
        enum isValidOperatorImpl = isValidOperatorImpl!(R, op, L);
    else
    {
        static if(hasStaticRows!L)
        {
            static if(hasStaticRows!R)
                enum _isValidR = L.rows == R.rows;
            else
                enum _isValidR = true;
        }
        else
            enum _isValidR = true;

        static if(hasStaticColumns!L)
        {
            static if(hasStaticColumns!R)
                enum _isValidC = L.cols == R.cols;
            else
                enum _isValidC = true;
        }
        else
            enum _isValidC = true;

        enum isValidOperatorImpl = _isValidR && _isValidC;
    }


    struct Inferred(M, size_t r, size_t c)
    if(r == wild || c == wild)
    {
        enum size_t rows = M.inferSize(r, c).rows;
        enum size_t cols = M.inferSize(r, c).cols;

        auto opIndex(size_t i, size_t j){ return M.init[i, j]; }
    }
}


template isValidOperatorImpl(L, string op, R)
if(isMatrix!L && isMatrix!R && op == "*")
{
    struct Inferred(M, size_t r, size_t c)
    if(r == wild || c == wild)
    {
        enum size_t rows = M.inferSize(r, c).rows;
        enum size_t cols = M.inferSize(r, c).cols;

        auto opIndex(size_t i, size_t j){ return M.init[i, j]; }
    }


    static if(isInferableMatrix!L && isInferableMatrix!R)
        enum isValidOperatorImpl = false;
    else static if(isInferableMatrix!L)
    {
        static if(hasStaticRows!R)
            enum isValidOperatorImpl = isValidOperatorImpl!(Inferred!(L, wild, R.rows), op, R);
        else
            enum isValidOperatorImpl = true;
    }
    else static if(isInferableMatrix!R)
    {
        static if(hasStaticColumns!L)
            enum isValidOperatorImpl = isValidOperatorImpl!(L, op, Inferred!(R, L.cols, wild));
        else
            enum isValidOperatorImpl = true;
    }
    else
    {
        static if(hasStaticColumns!L && hasStaticRows!R)
            enum isValidOperatorImpl = L.cols == R.rows;
        else
            enum isValidOperatorImpl = true;
    }
}


template isValidOperatorImpl(L, string op, R)
if((isMatrix!L && !isMatrix!R) || (isMatrix!R && !isMatrix!L))
{
    static if(op != "+" && op != "-" && op != "*" && op != "/")
        enum isValidOperatorImpl = false;
    else
        enum isValidOperatorImpl = true;
}


unittest{
    static struct S(T, size_t r, size_t c){enum rows = r; enum cols = c; T opIndex(size_t i, size_t j){return T.init;}}
    alias Static1x1 = S!(int, 1, 1);
    alias Static1x2 = S!(int, 1, 2);
    alias Static2x1 = S!(int, 2, 1);
    alias Static2x2 = S!(int, 2, 2);

    static struct D(T){size_t rows = 1, cols = 1; T opIndex(size_t i, size_t j){return T.init;}}
    alias Dynamic = D!(int);

    static struct I(T){
        enum rows = 0, cols = 0;
        T opIndex(size_t i, size_t j){return T.init;}
        static InferredResult inferSize(size_t rs, size_t cs){
            if(isEqOrEitherEq0(rs, cs) && (rs != 0 || cs != 0))
                return InferredResult(true, max(rs, cs), max(rs, cs));
            else
                return InferredResult(false, 0, 0);
        }
    }
    alias Inferable = I!int;
    static assert(Inferable.inferSize(1, 0).isValid);

    alias T = Inferable;
    static assert(T.rows == wild);
    static assert(T.cols == wild);


    static assert( isValidOperator!(Static1x1, "+", Static1x1));
    static assert(!isValidOperator!(Static1x1, "+", Static1x2));
    static assert( isValidOperator!(Static1x2, "+", Static1x2));
    static assert(!isValidOperator!(Static1x2, "+", Static1x1));

    static assert( isValidOperator!(Static1x1, "+", Dynamic));
    static assert( isValidOperator!(Static1x2, "+", Dynamic));
    static assert( isValidOperator!(Dynamic, "+", Static1x1));
    static assert( isValidOperator!(Dynamic, "+", Static1x2));

    static assert( isValidOperator!(Static1x1, "+", Inferable));
    static assert(!isValidOperator!(Static1x2, "+", Inferable));
    static assert( isValidOperator!(Inferable, "+", Static1x1));
    static assert(!isValidOperator!(Inferable, "+", Static1x2));

    static assert( isValidOperator!(Static1x1, "*", Static1x1));
    static assert( isValidOperator!(Static1x1, "*", Static1x2));
    static assert(!isValidOperator!(Static1x2, "*", Static1x2));
    static assert(!isValidOperator!(Static1x2, "*", Static1x1));

    static assert( isValidOperator!(Static1x1, "*", Dynamic));
    static assert( isValidOperator!(Static1x2, "*", Dynamic));
    static assert( isValidOperator!(Dynamic, "*", Static1x1));
    static assert( isValidOperator!(Dynamic, "*", Static1x2));

    static assert( isValidOperator!(Static1x1, "*", Inferable));
    static assert( isValidOperator!(Static1x2, "*", Inferable));
    static assert( isValidOperator!(Inferable, "*", Static1x1));
    static assert( isValidOperator!(Inferable, "*", Static1x2));
}


/**
Expression Template Operator Species : 式テンプレート演算子の種類
*/
enum ETOSpec : size_t
{
    none = 0,
    matrixAddMatrix = (1 << 0),
    matrixSubMatrix = (1 << 1),
    matrixMulMatrix = (1 << 2),
    matrixAddScalar = (1 << 3),
    scalarAddMatrix = (1 << 4),
    matrixSubScalar = (1 << 5),
    scalarSubMatrix = (1 << 6),
    matrixMulScalar = (1 << 7),
    scalarMulMatrix = (1 << 8),
    matrixDivScalar = (1 << 9),
    scalarDivMatrix = (1 << 10),
    addScalar = (1 << 11),
    subScalar = (1 << 12),
    mulScalar = (1 << 13),
    divScalar = (1 << 14),
    opEquals = (1 << 15),
    toString = (1 << 16),
    modifiedOperator = addScalar | subScalar | mulScalar | divScalar,
    all = (1 << 17) -1,
}


/**
式テンプレートでの演算子の種類を返します。

Example:
---
static assert(ETOperatorSpec!(Matrix2i, "+", Matrix2i) == ETOSpec.matrixAddMatrix);
static assert(ETOperatorSpec!(Matrix2i, "-", Matrix2i) == ETOSpec.matrixSubMatrix);
static assert(ETOperatorSpec!(Matrix2i, "*", Matrix2i) == ETOSpec.matrixMulMatrix);
static assert(ETOperatorSpec!(Matrix2i, "*", int) == ETOSpec.matrixMulScalar);
static assert(ETOperatorSpec!(int, "*", Matrix1x2i) == ETOSpec.scalarMulMatrix);
---
*/
template ETOperatorSpec(A, string op, B)
if(isValidOperator!(A, op, B))
{
    static if(isMatrix!A && isMatrix!B)
        enum ETOSpec ETOperatorSpec = op == "+" ? ETOSpec.matrixAddMatrix
                                                : (op == "-" ? ETOSpec.matrixSubMatrix
                                                             : ETOSpec.matrixMulMatrix);
    else static if(isScalar!A)
        enum ETOSpec ETOperatorSpec = op == "+" ? ETOSpec.scalarAddMatrix
                                                : (op == "-" ? ETOSpec.scalarSubMatrix
                                                             : (op == "*" ? ETOSpec.scalarMulMatrix
                                                                          : ETOSpec.scalarDivMatrix));
    else
        enum ETOSpec ETOperatorSpec = op == "+" ? ETOSpec.matrixAddScalar
                                                : (op == "-" ? ETOSpec.matrixSubScalar
                                                             : (op == "*" ? ETOSpec.matrixMulScalar
                                                                          : ETOSpec.matrixDivScalar));
}
unittest{
    static struct S(T, size_t r, size_t c){enum rows = r; enum cols = c; T opIndex(size_t i, size_t j){return T.init;}}
    alias Matrix2i = S!(int, 2, 2);

    static assert(ETOperatorSpec!(Matrix2i, "+", Matrix2i) == ETOSpec.matrixAddMatrix);
    static assert(ETOperatorSpec!(Matrix2i, "-", Matrix2i) == ETOSpec.matrixSubMatrix);
    static assert(ETOperatorSpec!(Matrix2i, "*", Matrix2i) == ETOSpec.matrixMulMatrix);
    static assert(ETOperatorSpec!(Matrix2i, "*", int) == ETOSpec.matrixMulScalar);
    static assert(ETOperatorSpec!(int, "*", Matrix2i) == ETOSpec.scalarMulMatrix);
}


/**

*/
struct MatrixExpression(Lhs, string s, Rhs)
if(isValidOperator!(Lhs, s, Rhs) && (isInferableMatrix!Lhs && isInferableMatrix!Rhs) || (isInferableMatrix!Lhs && !isMatrix!Rhs) || (!isMatrix!Lhs && isInferableMatrix!Rhs))
{
    enum rows = wild;
    enum cols = wild;
    enum etoSpec = ETOperatorSpec!(Lhs, s, Rhs);


    static InferredResult inferSize(size_t r, size_t c)
    {
        static if(isInferableMatrix!Lhs && isInferableMatrix!Rhs)
        {
            static assert(s != "*");
            auto rLhs = Lhs.inferSize(r, c);
            auto rRhs = Rhs.inferSize(r, c);

            bool b = rLhs.isValid && rRhs.isValid && rLhs.rows == rRhs.rows && rLhs.cols == rRhs.cols;
            return InferredResult(b, rLhs.rows, rLhs.cols);
        }
        else static if(isInferableMatrix!Lhs)
            return Lhs.inferSize(r, c);
        else
            return Rhs.inferSize(r, c);
    }


    auto ref opIndex(size_t i, size_t j)
    {
      static if(etoSpec == ETOSpec.matrixAddMatrix)
        return this.lhs[i, j] + this.rhs[i, j];
      else static if(etoSpec == ETOSpec.matrixSubMatrix)
        return this.lhs[i, j] - this.rhs[i, j];
      else static if(etoSpec == ETOSpec.matrixMulMatrix)
      {
        static assert(0);
        return typeof(this.lhs[0, 0] + this.rhs[0, 0]).init;
      }
      else
      {
        static if(isMatrix!Lhs)
            return mixin("this.lhs[i, j] " ~ s ~ " this.rhs");
        else
            return mixin("this.lhs " ~ s ~ " this.rhs[i, j]");
      }
    }


    auto ref opIndex(size_t i, size_t j) const 
    {
      static if(etoSpec == ETOSpec.matrixAddMatrix)
        return this.lhs[i, j] + this.rhs[i, j];
      else static if(etoSpec == ETOSpec.matrixSubMatrix)
        return this.lhs[i, j] - this.rhs[i, j];
      else static if(etoSpec == ETOSpec.matrixMulMatrix)
      {
        static assert(0);
        return typeof(this.lhs[0, 0] + this.rhs[0, 0]).init;
      }
      else
      {
        static if(isMatrix!Lhs)
            return mixin("this.lhs[i, j] " ~ s ~ " this.rhs");
        else
            return mixin("this.lhs " ~ s ~ " this.rhs[i, j]");
      }
    }


    auto ref opIndex(size_t i, size_t j) immutable
    {
      static if(etoSpec == ETOSpec.matrixAddMatrix)
        return this.lhs[i, j] + this.rhs[i, j];
      else static if(etoSpec == ETOSpec.matrixSubMatrix)
        return this.lhs[i, j] - this.rhs[i, j];
      else static if(etoSpec == ETOSpec.matrixMulMatrix)
      {
        static assert(0);
        return typeof(this.lhs[0, 0] + this.rhs[0, 0]).init;
      }
      else
      {
        static if(isMatrix!Lhs)
            return mixin("this.lhs[i, j] " ~ s ~ " this.rhs");
        else
            return mixin("this.lhs " ~ s ~ " this.rhs[i, j]");
      }
    }


    mixin(defaultExprOps!(true));

  private:
    Lhs lhs;
    Rhs rhs;
}


struct MatrixExpression(Lhs, string s, Rhs)
if(isValidOperator!(Lhs, s, Rhs) && !((isInferableMatrix!Lhs && isInferableMatrix!Rhs) || (isInferableMatrix!Lhs && !isMatrix!Rhs) || (!isMatrix!Lhs && isInferableMatrix!Rhs)))
{
    enum etoSpec = ETOperatorSpec!(Lhs, s, Rhs);


    static if(isMatrix!Lhs && isMatrix!Rhs)
    {
        static if(s == "*")
        {
            static if(hasStaticRows!Lhs)
            {
                enum rows = Lhs.rows;
                private enum staticrows = rows;
            }
            else static if(isInferableMatrix!Lhs && hasStaticRows!Rhs)
            {
                enum rows = Lhs.inferSize(wild, Rhs.rows).rows;
                private enum staticrows = rows;
            }
            else static if(hasDynamicRows!Lhs)
            {
                @property size_t rows() const { return this.lhs.rows; }
                private enum staticrows = wild;
            }
            else static if(isInferableMatrix!Lhs && hasDynamicRows!Rhs)
            {
                @property size_t rows() const { return this.lhs.inferSize(wild, rhs.rows).rows; }
                private enum staticrows = wild;
            }
            else
                static assert(0);


            static if(hasStaticColumns!Rhs)
            {
                enum cols = Rhs.cols;
                private enum staticcols = cols;
            }
            else static if(isInferableMatrix!Rhs && hasStaticColumns!Lhs)
            {
                enum cols = Rhs.inferSize(Lhs.cols, wild).cols;
                private enum staticcols = cols;
            }
            else static if(hasDynamicRows!Rhs)
            {
                @property size_t cols() const { return this.rhs.cols; }
                private enum staticcols = wild;
            }
            else static if(isInferableMatrix!Rhs && hasDynamicColumns!Lhs)
            {
                @property size_t cols() const { return this.rhs.inferSize(lhs.cols, wild).cols; }
                private enum staticcols = wild;
            }
            else
                static assert(0);
        }
        else
        {
            static if(hasStaticRows!Lhs)
            {
                enum rows = Lhs.rows;
                private enum staticrows = rows;
            }
            else static if(hasStaticRows!Rhs)
            {
                enum rows = Rhs.rows;
                private enum staticrows = rows;
            }
            else static if(hasDynamicRows!Lhs)
            {
                @property size_t rows() const { return this.lhs.rows; }
                private enum staticrows = wild;
            }
            else
            {
                @property size_t rows() const { return this.rhs.rows; }
                private enum staticrows = wild;
            }

            static if(hasStaticColumns!Lhs)
            {
                enum cols = Lhs.cols;
                private enum staticcols = cols;
            }
            else static if(hasStaticColumns!Rhs)
            {
                enum cols = Rhs.cols;
                private enum staticcols = cols;
            }
            else static if(hasDynamicColumns!Lhs)
            {
                @property size_t cols() const { return this.lhs.cols; }
                private enum staticcols = wild;
            }
            else
            {
                @property size_t cols() const { return this.rhs.cols; }
                private enum staticcols = wild;
            }
        }
    }
    else static if(isMatrix!Lhs)
    {
        static assert(!isInferableMatrix!Lhs);

        static if(hasStaticRows!Lhs)
        {
            enum rows = Lhs.rows;
            private enum staticrows = rows;
        }
        else
        {
            @property size_t rows() const { return this.lhs.rows; }
            private enum staticrows = wild;
        }

        static if(hasStaticColumns!Lhs)
        {
            enum cols = Lhs.cols;
            private enum staticcols = cols;
        }
        else
        {
            @property size_t cols() const { return this.lhs.cols; }
            private enum staticcols = wild;
        }
    }
    else
    {
        static assert(!isInferableMatrix!Rhs);

        static if(hasStaticRows!Rhs)
        {
            enum rows = Rhs.rows;
            private enum staticrsght = rows;
        }
        else
        {
            @property size_t rows() const { return this.rhs.rows; }
            private enum staticrsght = wild;
        }

        static if(hasStaticColumns!Rhs)
        {
            enum cols = Rhs.cols;
            private enum staticcsght = cols;
        }
        else
        {
            @property size_t cols() const { return this.rhs.cols; }
            private enum staticcsght = wild;
        }
    }


    auto ref opIndex(size_t i, size_t j)
    in{
        assert(i < this.rows);
        assert(j < this.cols);
    }
    body{
      static if(etoSpec == ETOSpec.matrixAddMatrix)
        return this.lhs[i, j] + this.rhs[i, j];
      else static if(etoSpec == ETOSpec.matrixSubMatrix)
        return this.lhs[i, j] - this.rhs[i, j];
      else static if(etoSpec == ETOSpec.matrixMulMatrix)
      {
        Unqual!(typeof(this.lhs[0, 0] * this.rhs[0, 0])) sum = 0;

        static if(hasStaticColumns!Lhs)
            immutable cnt = Lhs.cols;
        else static if(hasStaticRows!Rhs)
            immutable cnt = Rhs.rows;
        else static if(hasDynamicColumns!Lhs)
            immutable cnt = this.lhs.cols;
        else
            immutable cnt = this.rhs.rows;

        foreach(k; 0 .. cnt)
            sum += this.lhs[i, k] * this.rhs[k, j];
        return sum;
      }
      else
      {
        static if(isMatrix!Lhs)
            return mixin("this.lhs[i, j] " ~ s ~ " this.rhs");
        else
            return mixin("this.lhs " ~ s ~ " this.rhs[i, j]");
      }
    }


    auto ref opIndex(size_t i, size_t j) const 
    in{
        assert(i < this.rows);
        assert(j < this.cols);
    }
    body{
      static if(etoSpec == ETOSpec.matrixAddMatrix)
        return this.lhs[i, j] + this.rhs[i, j];
      else static if(etoSpec == ETOSpec.matrixSubMatrix)
        return this.lhs[i, j] - this.rhs[i, j];
      else static if(etoSpec == ETOSpec.matrixMulMatrix)
      {
        Unqual!(typeof(this.lhs[0, 0] * this.rhs[0, 0])) sum = 0;

        static if(hasStaticColumns!Lhs)
            immutable cnt = Lhs.cols;
        else static if(hasStaticRows!Rhs)
            immutable cnt = Rhs.rows;
        else static if(hasDynamicColumns!Lhs)
            immutable cnt = this.lhs.cols;
        else
            immutable cnt = this.rhs.rows;

        foreach(k; 0 .. cnt)
            sum += this.lhs[i, k] * this.rhs[k, j];
        return sum;
      }
      else
      {
        static if(isMatrix!Lhs)
            return mixin("this.lhs[i, j] " ~ s ~ " this.rhs");
        else
            return mixin("this.lhs " ~ s ~ " this.rhs[i, j]");
      }
    }


    auto ref opIndex(size_t i, size_t j) immutable
    in{
        assert(i < this.rows);
        assert(j < this.cols);
    }
    body{
      static if(etoSpec == ETOSpec.matrixAddMatrix)
        return this.lhs[i, j] + this.rhs[i, j];
      else static if(etoSpec == ETOSpec.matrixSubMatrix)
        return this.lhs[i, j] - this.rhs[i, j];
      else static if(etoSpec == ETOSpec.matrixMulMatrix)
      {
        Unqual!(typeof(this.lhs[0, 0] * this.rhs[0, 0])) sum = 0;

        static if(hasStaticColumns!Lhs)
            immutable cnt = Lhs.cols;
        else static if(hasStaticRows!Rhs)
            immutable cnt = Rhs.rows;
        else static if(hasDynamicColumns!Lhs)
            immutable cnt = this.lhs.cols;
        else
            immutable cnt = this.rhs.rows;

        foreach(k; 0 .. cnt)
            sum += this.lhs[i, k] * this.rhs[k, j];
        return sum;
      }
      else
      {
        static if(isMatrix!Lhs)
            return mixin("this.lhs[i, j] " ~ s ~ " this.rhs");
        else
            return mixin("this.lhs " ~ s ~ " this.rhs[i, j]");
      }
    }


    mixin(defaultExprOps!(false));

  private:
    Lhs lhs;
    Rhs rhs;
}

/**
行列に関する演算子の式
*/
auto matrixExpression(string s, A, B)(auto ref A a, auto ref B b)
if(isValidOperator!(A, s, B))
{
    return MatrixExpression!(A, s, B)(a, b);
    static assert(isMatrix!(MatrixExpression!(A, s, B)));
}


/**
specにいれた種類の演算子をオーバーロードします。

Example:
---
struct Identity(T){
    enum rows = wild;
    enum cols = wild;

    T opIndex(size_t i, size_t j){
        return i == j ? cast(T)0 : cast(T)1;
    }

    mixin ExpressionTemplateOperators!(ETOSpec.all, wild, wild);
}
---
*/

template ExpressionOperators(size_t spec, size_t rs, size_t cs)
{
    
    enum stringMixin = 
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        static if(is(typeof({enum _unused_ = this.rows;}))) static assert(rows != 0);
    } ~


    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        static if(is(typeof({enum _unused_ = this.cols;})))static assert(cols != 0);
    } ~


    ( rs == 1 ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        alias cols length;
        alias length opDollar;


        auto ref opIndex(size_t i)
        in{
            assert(i < this.cols);
        }
        body{
            return this[0, i];
        }

      static if(is(typeof({this[0] = this[0];})))
      {
        auto ref opIndexAssign(S)(S value, size_t i)
        if(is(typeof(this[0, i] = value)))
        in{
            assert(i < this.cols);
        }
        body{
            return this[0, i] = value;
        }


        auto ref opIndexOpAssign(string op, S)(S value, size_t i)
        if(is(typeof(mixin("this[0, i] " ~ op ~ "= value"))))
        in{
            assert(i < this.cols);
        }
        body{
            return mixin("this[0, i] " ~ op ~ "= value");
        }
      }
    } : ( cs == 1 ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        alias rows length;
        alias length opDollar;


        auto ref opIndex(size_t i)
        in{
            assert(i < this.rows);
        }
        body{
            return this[i, 0];
        }

      static if(is(typeof({this[0] = this[0];})))
      {
        auto ref opIndexAssign(S)(S value, size_t i)
        if(is(typeof(this[i, 0] = value)))
        in{
            assert(i < this.rows);
        }
        body{
            return this[0, i] = value;
        }


        auto ref opIndexOpAssign(string op, S)(S value, size_t i)
        if(is(typeof(mixin("this[i, 0] " ~ op ~ "= value"))))
        in{
            assert(i < this.rows);
        }
        body{
            return mixin("this[i, 0] " ~ op ~ "= value");
        }
      }
    } : ""
    )) ~


    (rs == 1 && cs == 1 ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        auto det() @property
        {
            return this[0, 0];
        }

        alias det this;


        auto ref opAssign(S)(S value)
        if(is(typeof(this[0, 0] = value)))
        {
            return this[0, 0] = value;
        }


        auto ref opOpAssign(S)(S value)
        if(is(typeof(mixin("this[0, 0] " ~ op ~ "= value"))))
        {
            return mixin("this[0, 0] " ~ op ~ "= value");
        }
    } : ""
    ) ~


    (spec & ETOSpec.opEquals ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        bool opEquals(Rhs)(auto ref const Rhs mat)
        if(isMatrix!Rhs)
        {
            static assert(isValidOperator!(Unqual!(typeof(this)), "+", Rhs));

            static if(isInferableMatrix!Rhs)
            {
                auto result = Rhs.inferSize(this.rows, this.cols);
                if(!result.isValid)
                    return false;
            }
            else
            {
                if(this.rows != mat.rows)
                    return false;

                if(this.cols != mat.cols)
                    return false;
            }

            foreach(i; 0 .. this.rows)
                foreach(j; 0 .. this.cols)
                    if(this[i, j] != mat[i, j])
                        return false;
            return true;
        }
    } : ""
    ) ~


    (spec & ETOSpec.toString ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        /*   //dmd bug : toStringをONにすると、メモリをバカ食いする事象*/
        @property
        void toString(scope void delegate(const(char)[]) sink, string formatString) @system
        {
            //sink(formatString);
            //formattedWrite(sink, formatString.replace("%r$", "%2$").replace("%c$", "%3$"), this.toRange!(Major.row), this.rows, this.cols);
            foreach(i; 0 .. this.rows)
                foreach(j; 0 .. this.cols)
                    formattedWrite(sink, "%s, ", this[i, j]);
        }
    } : ""
    ) ~


    (spec & ETOSpec.matrixAddMatrix ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        auto opBinary(string op : "+", Rhs)(auto ref Rhs mat)
        if(isMatrix!Rhs)
        in{
            static if(isInferableMatrix!Rhs)
                assert(mat.inferSize(this.rows, this.cols).isValid);
            else
            {
                assert(mat.rows == this.rows);
                assert(mat.cols == this.cols);
            }
        }
        body{
            static assert(isValidOperator!(typeof(this), op, Rhs));
            return matrixExpression!"+"(this, mat);
        }
    } : ""
    ) ~


    (spec & ETOSpec.matrixSubMatrix ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        auto opBinary(string op : "-", Rhs)(auto ref Rhs mat)
        if(isMatrix!Rhs)
        in{
            static if(isInferableMatrix!Rhs)
                assert(mat.inferSize(this.rows, this.cols).isValid);
            else
            {
                assert(mat.rows == this.rows);
                assert(mat.cols == this.cols);
            }
        }
        body{
            static assert(isValidOperator!(typeof(this), op, Rhs));
            return matrixExpression!"-"(this, mat);
        }
    } : ""
    ) ~


    (spec & ETOSpec.matrixMulMatrix ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        auto opBinary(string op : "*", Rhs)(auto ref Rhs mat)
        if(isMatrix!Rhs)
        in{
            static if(isInferableMatrix!Rhs)
                assert(mat.inferSize(this.cols, wild).isValid);
            else
                assert(mat.rows == this.cols);
        }
        body{
            static assert(isValidOperator!(typeof(this), op, Rhs));
            return matrixExpression!"*"(this, mat);
        }
    } : ""
    ) ~


    (spec & ETOSpec.matrixAddScalar ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        auto opBinary(string op : "+", S)(S s)
        if(isScalar!S)
        {
            static assert(isValidOperator!(typeof(this), op, S));
            return matrixExpression!"+"(this, s);
        }
    } : ""
    ) ~


    (spec & ETOSpec.scalarAddMatrix ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        auto opBinaryRight(string op : "+", S)(S s)
        if(isScalar!S)
        {
            static assert(isValidOperator!(S, op, typeof(this)));
            return matrixExpression!"+"(s, this);
        }
    } : ""
    ) ~


    (spec & ETOSpec.matrixSubScalar ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        auto opBinary(string op : "-", S)(S s)
        if(isScalar!S)
        {
            static assert(isValidOperator!(typeof(this), op, S));
            return matrixExpression!"-"(this, s);
        }
    } : ""
    ) ~


    (spec & ETOSpec.scalarSubMatrix ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        auto opBinaryRight(string op : "-", S)(S s)
        if(isScalar!S)
        {
            static assert(isValidOperator!(S, op, typeof(this)));
            return matrixExpression!"-"(s, this);
        }
    } : ""
    ) ~


    (spec & ETOSpec.matrixMulScalar ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        auto opBinary(string op : "*", S)(S s)
        if(isScalar!S)
        {
            static assert(isValidOperator!(typeof(this), op, S));
            return matrixExpression!"*"(this, s);
        }
    } : ""
    ) ~


    (spec & ETOSpec.scalarMulMatrix ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        auto opBinaryRight(string op : "*", S)(S s)
        if(isScalar!S)
        {
            static assert(isValidOperator!(S, op, typeof(this)));
            return matrixExpression!"*"(s, this);
        }
    } : ""
    ) ~


    (spec & ETOSpec.matrixDivScalar ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        auto opBinary(string op : "/", S)(S s)
        if(isScalar!S)
        {
            static assert(isValidOperator!(typeof(this), op, S));
            return matrixExpression!"/"(this, s);
        }
    } : ""
    ) ~


    (spec & ETOSpec.scalarDivMatrix ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        auto opBinaryRight(string op : "/", S)(S s)
        if(isScalar!S)
        {
            static assert(isValidOperator!(S, op, typeof(this)));
            return matrixExpression!"/"(s, this);
        }
    } : ""
    ) ~ 


    (spec & ETOSpec.addScalar ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        void opOpAssign(string op : "+", S)(S scalar)
        if(is(typeof(this[0, 0] += scalar)))
        {
            foreach(r; 0 .. rows)
                foreach(c; 0 .. cols)
                    this[r, c] +=  scalar;
        }
    } : ""
    ) ~


    (spec & ETOSpec.subScalar ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        void opOpAssign(string op : "-", S)(S scalar)
        if(is(typeof(this[0, 0] -= scalar)))
        {
            foreach(r; 0 .. rows)
                foreach(c; 0 .. cols)
                    this[r, c] -=  scalar;
        }
    } : ""
    ) ~


    (spec & ETOSpec.mulScalar ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        void opOpAssign(string op : "*", S)(S scalar)
        if(is(typeof(this[0, 0] *= scalar)))
        {
            foreach(r; 0 .. rows)
                foreach(c; 0 .. cols)
                    this[r, c] *=  scalar;
        }
    } : ""
    ) ~


    (spec & ETOSpec.divScalar ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        void opOpAssign(string op : "/", S)(S scalar)
        if(is(typeof(this[0, 0] /= scalar)))
        {
            foreach(r; 0 .. rows)
                foreach(c; 0 .. cols)
                    this[r, c] /=  scalar;
        }
    } : ""
    );


    mixin template templateMixin()
    {
        mixin(stringMixin);
    }
}


//inferable matrixのため
template ExpressionOperatorsInferable(size_t spec)
{
    enum stringMixin = 
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
    bool opEquals(Rhs)(auto ref const Rhs mat)
    if(isMatrix!Rhs && !is(Unqual!Rhs == typeof(this)) && !isInferableMatrix!(Rhs))
    {
        static assert(isValidOperator!(Unqual!(typeof(this)), "+", Rhs));

        if(!this.inferSize(mat.rows, mat.cols).isValid)
            return false;

        foreach(i; 0 .. mat.rows)
            foreach(j; 0 .. mat.cols)
                if(this[i, j] != mat[i, j])
                    return false;
        return true;
    }
    } ~


    (spec & ETOSpec.matrixAddMatrix ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        auto opBinary(string op : "+", Rhs)(auto ref Rhs mat)
        if(isMatrix!Rhs)
        in{
            static if(!isInferableMatrix!Rhs)
                assert(this.inferSize(mat.rows, mat.cols).isValid);
        }
        body{
            static assert(isValidOperator!(typeof(this), op, Rhs));
            return matrixExpression!"+"(this, mat);
        }
    } : ""
    ) ~


    (spec & ETOSpec.matrixSubMatrix ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        auto opBinary(string op : "-", Rhs)(auto ref Rhs mat)
        if(isMatrix!Rhs)
        in{
            static if(!isInferableMatrix!Rhs)
                assert(this.inferSize(mat.rows, mat.cols).isValid);
        }
        body{
            static assert(isValidOperator!(typeof(this), op, Rhs));
            return matrixExpression!"-"(this, mat);
        }
    } : ""
    ) ~


    (spec & ETOSpec.matrixMulMatrix ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        auto opBinary(string op : "*", Rhs)(auto ref Rhs mat)
        if(isMatrix!Rhs)
        in{
            static if(!isInferableMatrix!Rhs)
                assert(this.inferSize(wild, mat.rows).isValid);
        }
        body{
            static assert(isValidOperator!(typeof(this), op, Rhs));
            return matrixExpression!"*"(this, mat);
        }
    } : ""
    ) ~


    (spec & ETOSpec.matrixAddScalar ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        auto opBinary(string op : "+", S)(S s)
        if(isScalar!S)
        {
            static assert(isValidOperator!(typeof(this), op, S));
            return matrixExpression!"+"(this, s);
        }
    } : ""
    ) ~


    (spec & ETOSpec.scalarAddMatrix ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        auto opBinaryRight(string op : "+", S)(S s)
        if(isScalar!S)
        {
            static assert(isValidOperator!(S, op, typeof(this)));
            return matrixExpression!"+"(s, this);
        }
    } : ""
    ) ~


    (spec & ETOSpec.matrixSubScalar ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        auto opBinary(string op : "-", S)(S s)
        if(isScalar!S)
        {
            static assert(isValidOperator!(typeof(this), op, S));
            return matrixExpression!"-"(this, s);
        }
    } : ""
    ) ~


    (spec & ETOSpec.scalarSubMatrix ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        auto opBinaryRight(string op : "-", S)(S s)
        if(isScalar!S)
        {
            static assert(isValidOperator!(S, op, typeof(this)));
            return matrixExpression!"-"(s, this);
        }
    } : ""
    ) ~


    (spec & ETOSpec.matrixMulScalar ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        auto opBinary(string op : "*", S)(S s)
        if(isScalar!S)
        {
            static assert(isValidOperator!(typeof(this), op, S));
            return matrixExpression!"*"(this, s);
        }
    } : ""
    ) ~


    (spec & ETOSpec.scalarMulMatrix ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        auto opBinaryRight(string op : "*", S)(S s)
        if(isScalar!S)
        {
            static assert(isValidOperator!(S, op, typeof(this)));
            return matrixExpression!"*"(s, this);
        }
    } : ""
    ) ~


    (spec & ETOSpec.matrixDivScalar ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        auto opBinary(string op : "/", S)(S s)
        if(isScalar!S)
        {
            static assert(isValidOperator!(typeof(this), op, S));
            return matrixExpression!"/"(this, s);
        }
    } : ""
    ) ~

    (spec & ETOSpec.scalarDivMatrix ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        auto opBinaryRight(string op : "/", S)(S s)
        if(isScalar!S)
        {
            static assert(isValidOperator!(S, op, typeof(this)));
            return matrixExpression!"/"(s, this);
        }
    } : ""
    ) ~


    (spec & ETOSpec.addScalar ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        void opOpAssign(string op : "+", S)(S scalar)
        if(isScalar!S && is(typeof(this[0, 0] += scalar)))
        {
            foreach(r; 0 .. rows)
                foreach(c; 0 .. cols)
                    this[r, c] +=  scalar;
        }
    } : ""
    ) ~


    (spec & ETOSpec.subScalar ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        void opOpAssign(string op : "-", S)(S scalar)
        if(isScalar!S && is(typeof(this[0, 0] -= scalar)))
        {
            foreach(r; 0 .. rows)
                foreach(c; 0 .. cols)
                    this[r, c] -=  scalar;
        }
    } : ""
    ) ~


    (spec & ETOSpec.mulScalar ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        void opOpAssign(string op : "*", S)(S scalar)
        if(isScalar!S && is(typeof(this[0, 0] *= scalar)))
        {
            foreach(r; 0 .. rows)
                foreach(c; 0 .. cols)
                    this[r, c] *=  scalar;
        }
    } : ""
    ) ~


    (spec & ETOSpec.divScalar ?
    format("#line %s \"%s\"", __LINE__+1, __FILE__) ~ q{
        void opOpAssign(string op : "+", S)(S scalar)
        if(isScalar!S && is(typeof(this[0, 0] /= scalar)))
        {
            foreach(r; 0 .. rows)
                foreach(c; 0 .. cols)
                    this[r, c] /=  scalar;
        }
    } : ""
    );


    mixin template templateMixin()
    {
        mixin(stringMixin);
    }
}


/**

*/
template defaultExprOps(bool isInferable = false)
{
    enum defaultExprOps = 
    isInferable ?
    q{
        //mixin(ExpressionOperatorsInferable!(ETOSpec.all & ~ETOSpec.opEquals & ~ETOSpec.toString).stringMixin);
        mixin(ExpressionOperatorsInferable!(ETOSpec.modifiedOperator).stringMixin);
        const{mixin(ExpressionOperatorsInferable!(ETOSpec.all & ~ETOSpec.modifiedOperator).stringMixin);}
        //immutable{mixin(ExpressionOperatorsInferable!(ETOSpec.all & ~ETOSpec.opEquals & ~ETOSpec.toString).stringMixin);}
    }
    :
    q{
        //mixin(ExpressionOperators!(ETOSpec.all & ~ETOSpec.opEquals & ~ETOSpec.toString, mixin(is(typeof({enum _unused_ = rows;})) ? "this.rows" : "wild"), mixin(is(typeof({enum _unused_ = cols;})) ? "this.cols" : "wild")).stringMixin);
        mixin(ExpressionOperators!(ETOSpec.modifiedOperator, mixin(is(typeof({enum _unused_ = rows;})) ? "this.rows" : "wild"), mixin(is(typeof({enum _unused_ = cols;})) ? "this.cols" : "wild")).stringMixin);
        const{mixin(ExpressionOperators!(ETOSpec.all & ~ETOSpec.modifiedOperator, mixin(is(typeof({enum _unused_ = rows;})) ? "this.rows" : "wild"), mixin(is(typeof({enum _unused_ = cols;})) ? "this.cols" : "wild")).stringMixin);}
        //immutable{mixin(ExpressionOperators!(ETOSpec.all & ~ETOSpec.opEquals & ~ETOSpec.toString, mixin(is(typeof({enum _unused_ = rows;})) ? "this.rows" : "wild"), mixin(is(typeof({enum _unused_ = cols;})) ? "this.cols" : "wild")).stringMixin);}
    };
}


unittest{
    scope(failure) {writefln("Unittest failure :%s(%s)", __FILE__, __LINE__); stdout.flush();}
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}


    static struct M(size_t r, size_t c)
    {
        enum rows = r;
        enum cols = c;

        size_t opIndex(size_t i, size_t j) inout {return i + j;}

        //inout:
        mixin(defaultExprOps!(false));
    }

    alias S3 = M!(3, 3);
    alias S23 = M!(2, 3);

    static assert(isMatrix!S3);
    static assert(hasStaticRows!S3);
    static assert(hasStaticColumns!S3);
    static assert(isMatrix!S23);
    static assert(hasStaticRows!S23);
    static assert(hasStaticColumns!S23);


    M!(1, 1) m11;
    static assert(isMatrix!(typeof(m11)));
    int valueM11 = m11;

    M!(3, 1) m31;
    M!(1, 3) m13;
    valueM11 = m13 * m31;


    static struct I{
        enum rows = wild;
        enum cols = wild;

        size_t opIndex(size_t i, size_t j) inout { return i == j ? 1  : 0;}

        static InferredResult inferSize(size_t r, size_t c)
        {
            if(r == wild && c == wild)
                return InferredResult(false);
            else if(isEqOrEitherEq0(r, c))
                return InferredResult(true, max(r, c), max(r, c));
            else
                return InferredResult(false);
        }

        mixin(defaultExprOps!(true));
    }

    static assert(isMatrix!I);
    static assert(isInferableMatrix!I);
    static assert( I.inferSize(0, 1).isValid);
    static assert( I.inferSize(3, 3).isValid);
    static assert(!I.inferSize(1, 3).isValid);


    static struct D{
        size_t rows;
        size_t cols;

        size_t opIndex(size_t i, size_t j) inout {return i + j;}

        mixin(defaultExprOps!(false));
    }
    static assert(isMatrix!D);
    static assert(hasDynamicRows!D);
    static assert(hasDynamicColumns!D);


    S3 a;
    auto add = a + a;
    static assert(isMatrix!(typeof(add)));
    static assert(hasStaticRows!(typeof(add)));
    static assert(hasStaticColumns!(typeof(add)));
    assert(add[0, 0] == 0); assert(add[0, 1] == 2); assert(add[0, 2] == 4);
    assert(add[1, 0] == 2); assert(add[1, 1] == 4); assert(add[1, 2] == 6);
    assert(add[2, 0] == 4); assert(add[2, 1] == 6); assert(add[2, 2] == 8);

    auto mul = a * a;
    static assert(isMatrix!(typeof(mul)));
    static assert(hasStaticRows!(typeof(mul)));
    static assert(hasStaticColumns!(typeof(mul)));
    assert(mul[0, 0] == 5); assert(mul[0, 1] == 8); assert(mul[0, 2] ==11);
    assert(mul[1, 0] == 8); assert(mul[1, 1] ==14); assert(mul[1, 2] ==20);
    assert(mul[2, 0] ==11); assert(mul[2, 1] ==20); assert(mul[2, 2] ==29);

    auto sadd = a + 3;
    static assert(isMatrix!(typeof(sadd)));
    static assert(hasStaticRows!(typeof(sadd)));
    static assert(hasStaticColumns!(typeof(sadd)));
    assert(sadd[0, 0] == 3); assert(sadd[0, 1] == 4); assert(sadd[0, 2] == 5);
    assert(sadd[1, 0] == 4); assert(sadd[1, 1] == 5); assert(sadd[1, 2] == 6);
    assert(sadd[2, 0] == 5); assert(sadd[2, 1] == 6); assert(sadd[2, 2] == 7);

    // auto addMasS = a + m11;

    auto add5 = a + a + cast(const)(a) * 3;
    static assert(isMatrix!(typeof(add5)));
    static assert(hasStaticRows!(typeof(add5)));
    static assert(hasStaticColumns!(typeof(add5)));
    assert(add5 == a * 5);

    I i;
    auto addi = a + i;
    static assert(isMatrix!(typeof(addi)));
    static assert(hasStaticRows!(typeof(addi)));    static assert(typeof(addi).rows == 3);
    static assert(hasStaticColumns!(typeof(addi))); static assert(typeof(addi).cols == 3);
    assert(addi[0, 0] == 1); assert(addi[0, 1] == 1); assert(addi[0, 2] == 2);
    assert(addi[1, 0] == 1); assert(addi[1, 1] == 3); assert(addi[1, 2] == 3);
    assert(addi[2, 0] == 2); assert(addi[2, 1] == 3); assert(addi[2, 2] == 5);

    auto i2 = i * 2;
    static assert(isMatrix!(typeof(i2)));
    static assert(isInferableMatrix!(typeof(i2)));
    static assert( typeof(i2).inferSize(0, 1).isValid);
    static assert( typeof(i2).inferSize(3, 3).isValid);
    static assert(!typeof(i2).inferSize(1, 3).isValid);

    auto addi2 = a + i2;
    static assert(isMatrix!(typeof(addi2)));
    static assert(hasStaticRows!(typeof(addi2)));    static assert(typeof(addi2).rows == 3);
    static assert(hasStaticColumns!(typeof(addi2))); static assert(typeof(addi2).cols == 3);
    assert(addi2[0, 0] == 2); assert(addi2[0, 1] == 1); assert(addi2[0, 2] == 2);
    assert(addi2[1, 0] == 1); assert(addi2[1, 1] == 4); assert(addi2[1, 2] == 3);
    assert(addi2[2, 0] == 2); assert(addi2[2, 1] == 3); assert(addi2[2, 2] == 6);

    static assert(!is(typeof(S23.init + i)));
    static assert(!is(typeof(i + S23.init)));
    assert(S23.init * i == S23.init);
    assert(i * S23.init == S23.init);


    import core.exception, std.exception;

    D d33 = D(3, 3);
    auto addsd = a + d33;
    static assert(isMatrix!(typeof(addsd)));
    static assert(hasStaticRows!(typeof(addsd)));
    static assert(hasStaticColumns!(typeof(addsd)));
    assert(addsd == a * 2);
    assert(addsd == d33 * 2);

    auto addsdr = d33 + a;
    static assert(isMatrix!(typeof(addsdr)));
    static assert(hasStaticRows!(typeof(addsdr)));
    static assert(hasStaticColumns!(typeof(addsdr)));
    assert(addsdr == addsd);
    assert(addsdr == addsd);

    assert(collectException!AssertError(D(2, 3) + a));
    assert(collectException!AssertError(D(2, 3) + i));
    assert(D(2, 3) * i == D(2, 3));

    assert(collectException!AssertError(D(2, 3) + D(2, 2)));
    assert(collectException!AssertError(D(2, 3) + D(3, 3)));
    assert((D(2, 3) + D(2, 3)).rows == 2);
    assert((D(2, 3) + D(2, 3)).cols == 3);
    assert((D(2, 3) * D(3, 4)).rows == 2);
    assert((D(2, 3) * D(3, 4)).cols == 4);

    auto mulds = d33 * 3;
    assert(mulds == d33 + d33 + d33);
}


/**
InferableMatrixをある大きさの行列へ固定します。
*/
auto congeal(size_t r, size_t c, A)(auto ref A mat)
if(isInferableMatrix!A && A.inferSize(r, c).isValid)
{
    static struct Result()
    {
        enum size_t rows = A.inferSize(r, c).rows;
        enum size_t cols = A.inferSize(r, c).cols;


        auto ref opIndex(size_t i, size_t j) inout
        in{
            assert(i < rows);
            assert(j < cols);
        }
        body{
            return _mat[i, j];
        }


      static if(gigue.core.hasAssignableElements!A)
      {
        void opIndexAssign(E)(E e, size_t i, size_t j)
        in{
            assert(i < rows);
            assert(j < cols);
        }
        body{
            _mat[i, j] = e;
        }
      }

        mixin(defaultExprOps!(false));

      private:
        A _mat;
    }

    return Result!()(mat);
}
unittest{
    scope(failure) {writefln("Unittest failure :%s(%s)", __FILE__, __LINE__); stdout.flush();}
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}


    static struct I{
        enum rows = wild;
        enum cols = wild;

        size_t opIndex(size_t i, size_t j) inout { return i == j ? 1  : 0;}

        static InferredResult inferSize(size_t r, size_t c)
        {
            if(r == wild && c == wild)
                return InferredResult(false);
            else if(isEqOrEitherEq0(r, c))
                return InferredResult(true, max(r, c), max(r, c));
            else
                return InferredResult(false);
        }

        mixin(defaultExprOps!(true));
    }

    static assert(isMatrix!I);
    static assert(isInferableMatrix!I);
    static assert( I.inferSize(0, 1).isValid);
    static assert( I.inferSize(3, 3).isValid);
    static assert(!I.inferSize(1, 3).isValid);

    I id;
    auto i3x3 = id.congeal!(3, 3)();
    static assert(isMatrix!(typeof(i3x3)));
    static assert(hasStaticRows!(typeof(i3x3)));
    static assert(hasStaticColumns!(typeof(i3x3)));
    static assert(i3x3.rows == 3);
    static assert(i3x3.cols == 3);
    assert(i3x3[0, 0] == 1);assert(i3x3[1, 0] == 0);assert(i3x3[2, 0] == 0);
    assert(i3x3[0, 1] == 0);assert(i3x3[1, 1] == 1);assert(i3x3[2, 1] == 0);
    assert(i3x3[0, 2] == 0);assert(i3x3[1, 2] == 0);assert(i3x3[2, 2] == 1);
}


auto congeal(A)(auto ref A mat, size_t r, size_t c)
if(isInferableMatrix!A)
in{
    assert(A.inferSize(r, c).isValid);
}
body{
    static struct Result()
    {
        @property size_t rows() const { return _rs; }
        @property size_t cols() const { return _cs; }

        auto ref opIndex(size_t i, size_t j) inout
        in{
            assert(i < rows);
            assert(j < cols);
        }
        body{
            return _mat[i, j];
        }


      static if(gigue.core.hasAssignableElements!A)
      {
        void opIndexAssign(E)(E e, size_t i, size_t j)
        in{
            assert(i < rows);
            assert(j < cols);
        }
        body{
            _mat[i, j] = e;
        }
      }

        mixin(defaultExprOps!(false));

        size_t _rs, _cs;
        A _mat;
    }


    return Result!()(r, c, mat);
}
unittest{
    scope(failure) {writefln("Unittest failure :%s(%s)", __FILE__, __LINE__); stdout.flush();}
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}


    static struct I{
        enum rows = wild;
        enum cols = wild;

        size_t opIndex(size_t i, size_t j) inout { return i == j ? 1  : 0;}

        static InferredResult inferSize(size_t r, size_t c)
        {
            if(r == wild && c == wild)
                return InferredResult(false);
            else if(isEqOrEitherEq0(r, c))
                return InferredResult(true, max(r, c), max(r, c));
            else
                return InferredResult(false);
        }

        mixin(defaultExprOps!(true));
    }

    static assert(isMatrix!I);
    static assert(isInferableMatrix!I);
    static assert( I.inferSize(0, 1).isValid);
    static assert( I.inferSize(3, 3).isValid);
    static assert(!I.inferSize(1, 3).isValid);

    I id;
    auto i3x3 = id.congeal(3, 3);
    static assert(isMatrix!(typeof(i3x3)));
    static assert(hasDynamicRows!(typeof(i3x3)));
    static assert(hasDynamicColumns!(typeof(i3x3)));
    assert(i3x3.rows == 3);
    assert(i3x3.cols == 3);
    assert(i3x3[0, 0] == 1);assert(i3x3[1, 0] == 0);assert(i3x3[2, 0] == 0);
    assert(i3x3[0, 1] == 0);assert(i3x3[1, 1] == 1);assert(i3x3[2, 1] == 0);
    assert(i3x3[0, 2] == 0);assert(i3x3[1, 2] == 0);assert(i3x3[2, 2] == 1);
}


/**
DynamicMatrixをある大きさへ固定します。
*/
auto congeal(size_t r, size_t c, A)(auto ref A mat)
if(!isInferableMatrix!A && (hasDynamicRows!A || r == A.rows) || (hasDynamicColumns!A || c == A.cols))
in{
    assert(mat.rows == r);
    assert(mat.cols == c);
}
body{
    static struct Result()
    {
        enum rows = r;
        enum cols = c;


        auto ref opIndex(size_t i, size_t j) inout
        in{
            assert(i < rows);
            assert(j < cols);
        }
        body{
            return _mat[i, j];
        }


        static if(gigue.core.hasAssignableElements!A)
        {
          void opIndexAssign(E)(E e, size_t i, size_t j)
          in{
              assert(i < rows);
              assert(j < cols);
          }
          body{
              _mat[i, j] = e;
          }
        }


        mixin(defaultExprOps!(false));

      private:
        A _mat;
    }

    return Result!()(mat);
}
unittest{
    scope(failure) {writefln("Unittest failure :%s(%s)", __FILE__, __LINE__); stdout.flush();}
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}


    static struct D{
        size_t rows;
        size_t cols;

        size_t opIndex(size_t i, size_t j) inout {return i + j;}

        mixin(defaultExprOps!(false));
    }
    static assert(isMatrix!D);
    static assert(hasDynamicRows!D);
    static assert(hasDynamicColumns!D);

    D d3x3 = D(3, 3);
    auto s3x3 = d3x3.congeal!(3, 3)();
    static assert(isMatrix!(typeof(s3x3)));
    static assert(s3x3.rows == 3);
    static assert(s3x3.cols == 3);
    assert(s3x3[0, 0] == 0); assert(s3x3[1, 0] == 1); assert(s3x3[2, 0] == 2);
    assert(s3x3[0, 1] == 1); assert(s3x3[1, 1] == 2); assert(s3x3[2, 1] == 3);
    assert(s3x3[0, 2] == 2); assert(s3x3[1, 2] == 3); assert(s3x3[2, 2] == 4);
}


auto matrix(size_t r, size_t c, alias f, S...)(S sizes)
if(is(typeof(f(0, 0))) && S.length == ((r == 0) + (c == 0)) && (S.length == 0 || is(CommonType!S : size_t)))
{
    static struct Result()
    {
      static if(r)
        enum rows = r;
      else{
        private size_t _rs;
        @property size_t rows() const { return _rs; }
      }

      static if(c)
        enum cols = c;
      else{
        private size_t _cs;
        @property size_t cols() const { return _cs; }
      }


        auto ref opIndex(size_t i, size_t j) const
        in{
            assert(i < rows);
            assert(j < cols);
        }
        body{
            return f(i, j);
        }


        mixin(defaultExprOps!(false));
    }

    static if(S.length == 0)
      return Result!()();
    else static if(S.length == 1)
      return Result!()(sizes[0]);
    else static if(S.length == 2)
      return Result!()(sizes[0], sizes[1]);
}


auto matrix(alias f)(size_t r, size_t c)
if(is(typeof(f(0, 0))))
{
    return matrix!(0, 0, f)(r, c);
}


unittest{
    auto m = matrix!(2, 3, (i, j) => i * 3 + j);
    static assert(isMatrix!(typeof(m)));
    static assert(hasStaticRows!(typeof(m)));
    static assert(hasStaticColumns!(typeof(m)));
    static assert(m.rows == 2);
    static assert(m.cols == 3);
    assert(m[0, 0] == 0); assert(m[0, 1] == 1); assert(m[0, 2] == 2);
    assert(m[1, 0] == 3); assert(m[1, 1] == 4); assert(m[1, 2] == 5);


    auto md = matrix!((i, j) => i * 3 + j)(2, 3);
    assert(m == md);
}


auto matrix(alias f)()
{
    static struct Result()
    {
        enum rows = 0;
        enum cols = 0;


        static InferredResult inferSize(size_t r, size_t c)
        {
            if(r == 0 && c == 0)
                return  InferredResult(false);
            else if(r == 0 || c == 0)
                return InferredResult(true, max(r, c), max(r, c));
            else
                return InferredResult(true, r, c);
        }


        auto ref opIndex(size_t i, size_t j) const { return f(i, j); }


        mixin(defaultExprOps!(true));
    }

      return Result!()();
}
unittest{
    auto mi = matrix!((i, j) => i * 3 + j);
    static assert(isInferableMatrix!(typeof(mi)));
    assert(mi == matrix!((i, j) => i * 3 + j)(2, 3));
}


/**
内部で2次元配列を持つ行列です
*/
struct RefMatrix(T, size_t r = 0, size_t c = 0, Major mjr = Major.row)
{
    enum bool isColumnMajor = mjr == Major.column;
    enum bool isRowMajor    = mjr == Major.row;
    enum major = mjr;


    static if(r)
        enum size_t rows = r;
    else{
        @property
        size_t rows() const
        {
          static if(c == 1)
            return _array.length;
          else static if(isColumnMajor)
          {
            if(_array.length)
                return _array[0].length;
            else
                return 0;
          }
          else
            return _array.length;
        }


        @property
        void rows(size_t newSize)
        {
          static if(r == 1)
            return _array.length;
          else static if(isColummnMajor)
          {
            foreach(ref e; _array)
              e.length = newSize;
          }
          else
            _array.length = newSize;
        }
    }

    static if(c)
        enum size_t cols = c;
    else{
        @property
        size_t cols() const
        {
          static if(r == 1)
            return _array.length;
          else static if(isRowMajor){
            if(_array.length)
                return _array[0].length;
            else
                return 0;
          }
          else
            return _array.length;
        }

        @property
        void cols(size_t newSize)
        {
          static if(r == 1)
            return _array.length = newSize;
          else static if(isRowMajor)
          {
            foreach(ref e; _array)
              e.length = newSize;
          }
          else
            _array.length = newSize;
        }
    }


    auto ref opIndex(size_t i, size_t j) inout
    in{
        assert(i < rows);
        assert(j < cols);
    }
    body{
      static if(r == 1 || c == 1)
        return _array[i+j];
      else static if(isRowMajor)
        return _array[i][j];
      else
        return _array[j][i];
    }


    mixin(defaultExprOps!(false));


    void opAssign(M)(auto ref M mat)
    if(!is(typeof(this) == M) && isValidOperator!(typeof(this), "+", M))
    in{
      static if(r && hasDynamicRows!M)
        assert(this.rows == mat.rows);

      static if(c && hasDynamicColumns!M)
        assert(this.cols == mat.cols);
    }
    body
    {
        immutable rl = mat.rows,
                  cl = mat.cols,
                  n = rl * cl;

        if(_buffer.length < n)
            _buffer.length = n;

        foreach(i; 0 .. rl)
            foreach(j; 0 .. cl)
            {
                static if(major == Major.row)
                    _buffer[i * cl + j] = mat[i, j];
                else
                    _buffer[j * rl + i] = mat[i, j];
            }

      static if(!r)
        this.rows = rl;

      static if(!c)
        this.cols = cl;

      static if(r == 1 || c == 1)
        _array[] = _buffer[0 .. n][];
      else{
        foreach(i; 0 .. _array.length)
            _array[i][] = _buffer[i * _array[0].length .. (i+1) * _array[0].length][];
      }
    }


    @property
    void noAlias(M)(auto ref M mat)
    if(isValidOperator!(typeof(this), "+", M) && ((!r || !hasStaticRows!M) || is(typeof({static assert(this.rows == M.rows);})))
                                              && ((!c || !hasStaticColumns!M) || is(typeof({static assert(this.cols == M.cols);}))))
    in{
      static if(r)
        assert(this.rows == mat.rows);

      static if(c)
        assert(this.cols == mat.cols);
    }
    body{
      static if(is(typeof(this) == M))
        this = mat;
      else
      {
        immutable rl = mat.rows,
                  cl = mat.cols,
                  n = rl * cl;

        if(_array.length < n)
            _array.length = n;

        static if(!r)
          this.rows = rl;

        static if(!c)
          this.cols = cl;

        foreach(i; 0 .. rl)
            foreach(j; 0 .. cl)
            {
              static if(r == 1 || c == 1)
                _array[i+j] = mat[i, j];
              static if(major == Major.row)
                  _array[i][j] = mat[i, j];
              else
                  _array[j][i] = mat[i, j];
            }
      }
    }


    @property
    auto reference() inout
    {
        return this;
    }


  private:
  static if(r == 1 || c == 1)
    T[] _array;
  else
    T[][] _array;

  static:
    T[] _buffer;
}


RefMatrix!(T, r, c, mjr) matrix(size_t r, size_t c, Major mjr = Major.row, T)(T[] arr)
if(r != 0 || c != 0)
in{
  static if(r != 0 && c != 0)
    assert(arr.length == r * c);
  else{
    assert(!(arr.length % r + c));
  }
}
body{
    static if(r == 1 || c == 1)
      return typeof(return)(arr);
    else
    {
      immutable rs = r == 0 ? arr.length / c : r;
      immutable cs = c == 0 ? arr.length / r : c;

      immutable h = mjr == Major.row ? rs : cs,
                w = mjr == Major.row ? cs : rs;

      typeof(return) dst;
      dst._array.length = h;
      foreach(i; 0 .. h)
        dst._array[i] = arr[i * w .. (i+1) * w];

      return dst;
    }
}

unittest{
    scope(failure) {writefln("Unittest failure :%s(%s)", __FILE__, __LINE__); stdout.flush();}
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}

    auto mr = matrix!(2, 3)([0, 1, 2, 3, 4, 5]);
    static assert(isMatrix!(typeof(mr)));
    static assert(hasStaticRows!(typeof(mr)));
    static assert(hasStaticColumns!(typeof(mr)));
    static assert(mr.rows == 2);
    static assert(mr.cols == 3);
    assert(mr[0, 0] == 0); assert(mr[0, 1] == 1); assert(mr[0, 2] == 2);
    assert(mr[1, 0] == 3); assert(mr[1, 1] == 4); assert(mr[1, 2] == 5);

    mr[0, 0] = 10;
    assert(mr[0, 0] == 10); assert(mr[0, 1] == 1); assert(mr[0, 2] == 2);
    assert(mr[1, 0] ==  3); assert(mr[1, 1] == 4); assert(mr[1, 2] == 5);

    mr = matrix!(2, 3, (i, j) => (i*3+j)*2);
    assert(mr[0, 0] == 0); assert(mr[0, 1] == 2); assert(mr[0, 2] == 4);
    assert(mr[1, 0] == 6); assert(mr[1, 1] == 8); assert(mr[1, 2] == 10);

    auto mc = matrix!(2, 3, Major.column)([0, 1, 2, 3, 4, 5]);
    static assert(isMatrix!(typeof(mc)));
    static assert(hasStaticRows!(typeof(mc)));
    static assert(hasStaticColumns!(typeof(mc)));
    static assert(mc.rows == 2);
    static assert(mc.cols == 3);
    assert(mc[0, 0] == 0); assert(mc[0, 1] == 2); assert(mc[0, 2] == 4);
    assert(mc[1, 0] == 1); assert(mc[1, 1] == 3); assert(mc[1, 2] == 5);

    mc[0, 0] = 10;
    assert(mc[0, 0] == 10); assert(mc[0, 1] == 2); assert(mc[0, 2] == 4);
    assert(mc[1, 0] ==  1); assert(mc[1, 1] == 3); assert(mc[1, 2] == 5);

    mc = matrix!(2, 3, (i, j) => (i*3+j)*2);
    assert(mc[0, 0] == 0); assert(mc[0, 1] == 2); assert(mc[0, 2] == 4);
    assert(mc[1, 0] == 6); assert(mc[1, 1] == 8); assert(mc[1, 2] == 10);
}


RefMatrix!(T, wild, wild, mjr) matrix(Major mjr = Major.row, T)(T[] arr, size_t r, size_t c)
in{
  if(r != 0 && c != 0)
    assert(arr.length == r * c);
  else if(r != 0 || c != 0)
    assert(!(arr.length % (r + c)));
  else
    assert(0);
}
body{
    immutable rs = r == 0 ? arr.length / c : r;
    immutable cs = c == 0 ? arr.length / r : c;

    immutable h = mjr == Major.row ? rs : cs,
              w = mjr == Major.row ? rs : cs;

    typeof(return) dst;
    dst._array.length = h;
    foreach(e; 0 .. h)
      dst._array[i] = arr[i * w .. (i+1) * w];

    return dst;
}

unittest{
    scope(failure) {writefln("Unittest failure :%s(%s)", __FILE__, __LINE__); stdout.flush();}
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}

    auto mr = matrix!(2, 3)([0, 1, 2, 3, 4, 5]);
    static assert(isMatrix!(typeof(mr)));
    static assert(hasStaticRows!(typeof(mr)));
    static assert(hasStaticColumns!(typeof(mr)));
    static assert(mr.rows == 2);
    static assert(mr.cols == 3);
    assert(mr[0, 0] == 0); assert(mr[0, 1] == 1); assert(mr[0, 2] == 2);
    assert(mr[1, 0] == 3); assert(mr[1, 1] == 4); assert(mr[1, 2] == 5);


    auto mc = matrix!(2, 3, Major.column)([0, 1, 2, 3, 4, 5]);
    static assert(isMatrix!(typeof(mc)));
    static assert(hasStaticRows!(typeof(mc)));
    static assert(hasStaticColumns!(typeof(mc)));
    static assert(mc.rows == 2);
    static assert(mc.cols == 3);
    assert(mc[0, 0] == 0); assert(mc[0, 1] == 2); assert(mc[0, 2] == 4);
    assert(mc[1, 0] == 1); assert(mc[1, 1] == 3); assert(mc[1, 2] == 5);
}


auto matrix(Major mjr = Major.row, T, size_t N, size_t M)(auto ref T[M][N] arr)
{
  static if(mjr == Major.row)
    RefMatrix!(T, N, M, mjr) dst;
  else
    RefMatrix!(T, M, N, mjr) dst;

    dst._array.length = N;
    foreach(i; 0 .. N)
      dst._array[i] = arr[i];
    return dst;
}

unittest{
    scope(failure) {writefln("Unittest failure :%s(%s)", __FILE__, __LINE__); stdout.flush();}
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}

    int[3][2] arr = [[0, 1, 2], [3, 4, 5]];

    auto mr = matrix(arr);
    static assert(isMatrix!(typeof(mr)));
    static assert(hasStaticRows!(typeof(mr)));
    static assert(hasStaticColumns!(typeof(mr)));
    static assert(mr.rows == 2);
    static assert(mr.cols == 3);
    assert(mr[0, 0] == 0); assert(mr[0, 1] == 1); assert(mr[0, 2] == 2);
    assert(mr[1, 0] == 3); assert(mr[1, 1] == 4); assert(mr[1, 2] == 5);
}

auto matrix(Major mjr = Major.row, A)(A mat)
if(isMatrix!A && isInferableMatrix!A)
{
  static if(hasStaticRows!A && hasStaticColumns!A)
    RefMatrix!(ElementType!A, A.rows, A.cols, mjr) dst;
  else static if(hasStaticRows!A)
    RefMatrix!(ElementType!A, A.rows, wild, mjr) dst;
  else static if(hasStaticColumns!A)
    RefMatrix!(ElementType!A, wild, A.cols, mjr) dst;
  else
    RefMatrix!(ElementType!A, wild, wild, mjr) dst;

    dst.noAlias = mat;
    return dst;
}


auto matrix(size_t r, size_t c = 0, Major mjr = Major.row, A)(A mat)
if(isInferableMatrix!A && A.inferSize(r, c).isValid)
{
    return mat.congeal!(r, c).matrix();
}


auto matrix(Major mjr = Major.row, A)(A mat, size_t r, size_t c = 0)
if(isInferableMatrix!A)
in{
  assert(A.inferSize(r, c).isValid);
}
body{
  return mat.congeal(r, c).matrix();
}

unittest{

}


/**
要素がメモリ上に連続して存在するような行列
*/
struct Matrix(T, size_t r = 0, size_t c = 0, Major mjr = Major.row)
if(r != 0 && c != 0)
{
    enum bool isColumnMajor = mjr == Major.column;
    enum bool isRowMajor    = mjr == Major.row;
    enum major = mjr;
    enum size_t rows = r;
    enum size_t cols = c;


    this(M)(auto ref M mat)
    if(!is(M == typeof(this)))
    {
        this.noAlias = mat;
    }


    auto ref opIndex(size_t i, size_t j) inout
    in{
        assert(i < rows);
        assert(j < cols);
    }
    body{
        static if(major == Major.row)
            return _array[i * cols + j];
        else
            return _array[j * rows + i];
    }


    @property
    auto array() pure nothrow @safe inout
    {
        return _array[];
    }


    void opAssign(M)(auto ref M mat)
    if(isValidOperator!(typeof(this), "+", M))
    {
        foreach(i; 0 .. r)
            foreach(j; 0 .. c)
            {
                static if(major == Major.row)
                    Matrix._buffer[i * cols + j] = mat[i, j];
                else
                    Matrix._buffer[j * rows + i] = mat[i, j];
            }
        
        _array[] = Matrix._buffer[];
    }


    @property
    void noAlias(M)(auto ref M mat)
    if(isValidOperator!(typeof(this), "+", M))
    {
        foreach(i; 0 .. r)
            foreach(j; 0 .. c)
            {
                static if(major == Major.row)
                    _array[i * cols + j] = mat[i, j];
                else
                    _array[j * rows + i] = mat[i, j];
            }
    }


    @property
    auto reference()
    {
        return matrix!(r, c)(_array[]);
    }


    mixin(defaultExprOps!(false));


  private:
    T[r * c] _array;

  static:
    T[r * c] _buffer;
}

unittest{
    scope(failure) {writefln("Unittest failure :%s(%s)", __FILE__, __LINE__); stdout.flush();}
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}


    Matrix!(int, 3, 3) m;
    m[0, 0] = 0; m[0, 1] = 1; m[0, 2] = 2;
    m[1, 0] = 1; m[1, 1] = 2; m[1, 2] = 3;
    m[2, 0] = 2; m[2, 1] = 3; m[2, 2] = 4;

    Matrix!(int, 3, 3) m2 = m * m;
    m = m * m;
    assert(m == m2);
}

unittest{
    scope(failure) {writefln("Unittest failure :%s(%s)", __FILE__, __LINE__); stdout.flush();}
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}


    Matrix!(int, 2, 2, Major.row) mr;    // 2x2, int型, 行優先
    assert(mr.array.equal([0, 0, 0, 0]));       // 初期化される

    mr[0, 1] = 1;
    mr[1, 0] = 2;
    mr[1, 1] = 3;
    assert(mr.array.equal([0, 1, 2, 3]));       // 行優先順


    Matrix!(int, 2, 2, Major.column) mc; // 2x2, int型, 列優先
    assert(mc.array.equal([0, 0, 0, 0]));       // 初期化される

    mc[0, 1] = 1;
    mc[1, 0] = 2;
    mc[1, 1] = 3;
    assert(mc.array.equal([0, 2, 1, 3]));       // 列優先順


    Matrix!(int, 2, 2, Major.row) minit = mc;
    assert(minit.array.equal([0, 1, 2, 3]));   // 全要素12で初期化されている
}

unittest{
    scope(failure) {writefln("Unittest failure :%s(%s)", __FILE__, __LINE__); stdout.flush();}
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}


    Matrix!(int, 1, 3) m;
    m[0] = 3;
    assert(m[0] == 3);
    static assert(m.length == 3);
    assert(m[$-1] == 0);
}

unittest{
    scope(failure) {writefln("Unittest failure :%s(%s)", __FILE__, __LINE__); stdout.flush();}
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}


    Matrix!(int, 2, 2) m;
    auto rm = m.reference;

    assert(rm + rm == m + m);
    assert(rm - rm == m - m);
    assert(rm * rm == m * m);

    m[0, 0] = 1; m[0, 1] = 2;
    m[1, 0] = 2; m[1, 1] = 3;

    Matrix!(int, 2, 2) m2 = m;
    m.noAlias = m2 + m2;
    assert(m[0, 0] == 2);
}

unittest{
    scope(failure) {writefln("Unittest failure :%s(%s)", __FILE__, __LINE__); stdout.flush();}
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}


    alias Rvector!(int, 3) R;
    R rv1;
    static assert(rv1.rows == 1);
    static assert(rv1.cols == 3);
    static assert(rv1.length  == 3);

    rv1[0] = 3;
    assert(rv1[0] == 3);
    assert(rv1[1] == 0);
    assert(rv1[2] == 0);
    rv1.opOpAssign!"+"(3);
    assert(rv1[0] == 6);
    assert(rv1[1] == 3);
    assert(rv1[2] == 3);
    rv1[0] += 3;
    assert(rv1[0] == 9);
    assert(rv1[1] == 3);
    assert(rv1[2] == 3);

    Rvector!(int, 4) rv2;
    static assert(rv2.rows == 1);
    static assert(rv2.cols == 4);
    static assert(rv2.length  == 4);

    Cvector!(int, 3) cv1;
    static assert(cv1.rows == 3);
    static assert(cv1.cols == 1);
    static assert(cv1.length  == 3);

    Cvector!(int, 4) cv2;
    static assert(cv2.rows == 4);
    static assert(cv2.cols == 1);
    static assert(cv2.length  == 4);
}


///ditto
template Rvector(T, size_t size)
{
    alias Matrix!(T, 1, size, Major.row) Rvector;
}

///ditto
template Cvector(T, size_t size)
{
    alias Matrix!(T, size, 1, Major.column) Cvector;
}

///ditto
template Vector(T, size_t size)
{
    alias Cvector!(T, size) Vector;
}


/**
転置行列を返します。
*/
@property
auto transpose(A)(A mat)
if(isMatrix!A)
{
    static struct Transposed()
    {
      static if(isInferableMatrix!A)
      {
        enum size_t rows = wild;
        enum size_t cols = wild;

        static InferredResult inferSize(size_t r, size_t c)
        {
            return A.inferSize(c, r);
        }
      }
      else
      {
        static if(hasStaticColumns!A)
            enum size_t rows = A.cols;
        else
            @property auto ref rows() inout { return _mat.cols; }

        static if(hasStaticRows!A)
            enum size_t cols = A.rows;
        else
            @property auto ref cols() inout { return _mat.rows; }
      }


        auto ref opIndex(size_t i, size_t j) inout
        in{
            assert(i < rows);
            assert(j < cols);
        }
        body{
            return _mat[j, i];
        }


      static if(gigue.core.hasAssignableElements!A)
      {
        void opIndexAssign(E)(E e, size_t i, size_t j)
        in{
            assert(i < rows);
            assert(j < cols);
        }
        body{
            _mat[j, i] = e;
        }
      }

        mixin(defaultExprOps!(isInferableMatrix!A));


        auto ref transpose() @property inout
        {
            return _mat;
        }


      private:
        A _mat;
    }

    return Transposed!()(mat);
}
unittest{
    scope(failure) {writefln("Unittest failure :%s(%s)", __FILE__, __LINE__); stdout.flush();}
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}


    Matrix!(int, 2, 2) m;
    m[0, 0] = 0; m[0, 1] = 1;
    m[1, 0] = 2; m[1, 1] = 3;

    auto t = m.transpose;
    assert(t[0, 0] == 0);
    assert(t[0, 1] == 2);
    assert(t[1, 0] == 1);
    assert(t[1, 1] == 3);
}
unittest{
    scope(failure) {writefln("Unittest failure :%s(%s)", __FILE__, __LINE__); stdout.flush();}
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}


    Cvector!(int, 3) v;
    v[0] = 1; v[1] = 2; v[2] = 3;

    auto t = v.transpose;
    assert(t[0] == 1);
    assert(t[1] == 2);
    assert(t[2] == 3);

    assert(t.rows == 1);
    assert(t.cols == 3);


    static assert(is(typeof(v) == typeof(t.transpose)));
}

/**

*/
auto toRange(A)(A mat)
if(isMatrix!A)
{
    static struct ToRange()
    {
        static struct Element()
        {
            @property auto ref front() { return _mat[_r, _cf]; }

            @property auto ref back() { return _mat[_r, _cb]; }

            auto ref opIndex(size_t i) { i += _cf; return _mat[_r, i]; }

            void popFront() { ++_cf; }
            void popBack() { --_cb; }

            @property bool empty() { return _cf == _cb; }
            @property size_t length() { return _cb - _cf; }

            @property auto save() { return this; }

            auto opSlice() { return this.save; }
            auto opSlice(size_t i, size_t j) { return typeof(this)(_mat, _r, _cf + i, j); }


          private:
            A _mat;
            size_t _r;
            size_t _cf = 0;
            size_t _cb = A.cols;
        }


        @property auto front() { return Element!()(this._mat, _rf); }

        @property auto back() { return Element!()(this._mat, _rb); }

        auto opIndex(size_t i) { i += _rf; return Element!()(this._mat, i);}

        void popFront() { ++_rf; }
        void popBack() { --_rb; }

        @property bool empty() { return _rf == _rb; }
        @property size_t length() { return _rb - _rf; }

        @property auto save() { return this; }

        auto opSlice() { return this.save; }
        auto opSlice(size_t i, size_t j) { return typeof(this)(_mat, _rf + i, j); }

      private:
        A _mat;
        size_t _rf = 0;
        size_t _rb = A.rows;
    }


    return ToRange!()(mat);
}

unittest{
    scope(failure) {writefln("Unittest failure :%s(%s)", __FILE__, __LINE__); stdout.flush();}
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}


    Matrix!(int, 3, 3) rm33;
    rm33[0, 0] = 1; rm33[0, 1] = 2; rm33[0, 2] = 3;
    assert(equal!"equal(a, b)"(rm33.toRange, [[1, 2, 3], [0, 0, 0], [0, 0, 0]]));

    Matrix!(int, 1, 1) rm11;
    assert(equal!"equal(a, b)"(rm11.toRange, [[0]]));
}


/**
行列をレンジにします
*/
auto toFlatten(A)(A mat)
if(isMatrix!A && !isInferableMatrix!A)
{
    alias ElementType!A E;

    static struct ToFlatten()
    {
        @property
        auto ref front()
        {
            return _mat[_f / _mat.cols, _f % _mat.cols];
        }


        @property
        auto ref back()
        {
            return _mat[_b / _mat.cols, _b % _mat.cols];
        }


        auto ref opIndex(size_t i)
        in{
            assert(_f + i < _b);
        }body{
            i += _f;
            return _mat[i / _mat.cols, i % _mat.cols];
        }


        static if(hasAssignableElements!A)
        {
            @property
            void front(E v)
            {
                _mat[_f / _mat.cols, _f % _mat.cols] = v;
            }


            @property
            void back(E v)
            {
                _mat[_b / _mat.cols, _b % _mat.cols] = v;
            }


            void opIndexAssign(E v, size_t i)
            in{
                assert(_f + i < _b);
            }
            body{
                i += _f;
                _mat[i / _mat.cols, i % _mat.cols] = v;
            }
        }


        @property
        bool empty() pure nothrow @safe const
        {
            return _f >= _b;
        }


        void popFront() pure nothrow @safe
        {
            _f++;
        }


        void popBack() pure nothrow @safe
        {
            _b--;
        }


        @property
        size_t length() pure nothrow @safe const
        {
            return _b - _f;
        }


        alias length opDollar;


        @property
        typeof(this) save() pure nothrow @safe
        {
            return this;
        }


    private:
        A _mat;
        size_t _f, _b;
    }

    return ToFlatten!()(mat, 0, mat.rows * mat.cols);
}
unittest{
    scope(failure) {writefln("Unittest failure :%s(%s)", __FILE__, __LINE__); stdout.flush();}
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}


    Matrix!(int, 3, 3, Major.row) rm33;
    rm33[0, 0] = 1; rm33[0, 1] = 2; rm33[0, 2] = 3;
    //writeln(rm33.array.length;)
    assert(rm33.array == [1, 2, 3, 0, 0, 0, 0, 0, 0]);

    alias Rt1 = typeof(toFlatten(rm33));
    static assert(isRandomAccessRange!(Rt1));
    static assert(std.range.hasLvalueElements!(Rt1));
    static assert(std.range.hasAssignableElements!(Rt1));
    static assert(hasLength!(Rt1));
    assert(equal(rm33.toFlatten, rm33.array));

    Matrix!(int, 3, 3, Major.column) cm33;
    cm33[0, 0] = 1; cm33[0, 1] = 2; cm33[0, 2] = 3;
    assert(cm33.array == [1, 0, 0, 2, 0, 0, 3, 0, 0]);
    assert(equal(cm33.toFlatten, [1, 2, 3, 0, 0, 0, 0, 0, 0]));
}



/**
レンジから行列を作ります。
*/
auto toMatrix(size_t r, size_t c, Major mjr = Major.row, R)(R range)
if(isRandomAccessRange!R && isScalar!(Unqual!(std.range.ElementType!R)) && (mjr == Major.row ? (c != wild) : (r != wild)))
{
  //static if(mjr == Major.column)
    //return range.toMatrix!(c, r, Major.row).transpose;
  //else
  //{
    alias E = Unqual!(std.range.ElementType!R);

    static struct ToMatrix()
    {
      static if(r != wild)
        enum size_t rows = r;
      else
        auto ref rows() const @property
        {
            return _range.length / typeof(this).cols;
        }

      static if(c!= wild)
        enum size_t cols = c;
      else
        auto ref cols() const @property
        {
            return _range.length / typeof(this).rows;
        }


        auto ref opIndex(size_t i, size_t j) inout
        in{
            assert(i < rows || rows == 0);
            assert(j < cols || cols == 0);

          static if(hasLength!R && mjr == Major.row)
            assert(i * cols + j < this._range.length);

          static if(hasLength!R && mjr == Major.column)
            assert(j * rows + i < this._range.length);
        }
        body{
          static if(mjr == Major.row)
            return this._range[i * cols + j];
          else
            return this._range[j * rows + i];
        }


        mixin(defaultExprOps!(false));

      private:
        R _range;
    }

    return ToMatrix!()(range);
  //}
}

unittest{
    scope(failure) {writefln("Unittest failure :%s(%s)", __FILE__, __LINE__); stdout.flush();}
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}


    auto r = iota(4);
    auto mr = toMatrix!(2, 2)(r);
    assert(mr.toFlatten.equal([0, 1, 2, 3]));

    auto mc = r.toMatrix!(2, 2, Major.column);
    assert(mc.toFlatten.equal([0, 2, 1, 3]));
}
unittest{
    scope(failure) {writefln("Unittest failure :%s(%s)", __FILE__, __LINE__); stdout.flush();}
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}


    auto mem = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
    auto mr = mem.toMatrix!(3, 3);
    static assert(isMatrix!(typeof(mr)));
    static assert(hasLvalueElements!(typeof(mr)));
    assert(mr.toFlatten.equal(mem[0 .. 9]));

    mem.length = 16;
    auto mc = mem.toMatrix!(4, 4, Major.column);
    static assert(isMatrix!(typeof(mr)));
    static assert(hasLvalueElements!(typeof(mr)));
    mc[3, 3] = 15;
    assert(mc.transpose.toFlatten.equal([0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                                            0, 0, 0, 0, 0, 15]));
}
unittest{
    scope(failure) {writefln("Unittest failure :%s(%s)", __FILE__, __LINE__); stdout.flush();}
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}


    auto mem = [0, 1, 2, 3];
    auto r1 = mem.toMatrix!(wild, 1, Major.row);
    assert(equal(r1.toFlatten, mem[0 .. 4]));

    mem ~= [4, 5];
    auto r2 = mem.toMatrix!(wild, 2, Major.row);
    assert(r2[2, 0] == 4);
    assert(r2[2, 1] == 5);

    auto c1 = mem.toMatrix!(1, wild, Major.column);
    assert(equal(c1.toFlatten, mem));
}



auto toMatrix(size_t r, size_t c, Major mjr = Major.row, R)(R range)
if(isRandomAccessRange!R && isRandomAccessRange!(Unqual!(std.range.ElementType!R)) && isScalar!(Unqual!(std.range.ElementType!(Unqual!(std.range.ElementType!R)))))
{
    static struct ToMatrix()
    {
      static if(r != wild)
        enum size_t rows = r;
      else
        auto ref rows() const @property
        {
          static if(mjr == Major.row)
            return _range.length;
          else
            return _range[0].length;
        }


      static if(c != wild)
        enum size_t cols = c;
      else
        auto ref cols() const @property
        {
          static if(mjr == Major.column)
            return _range.length;
          else
            return _range[0].length;
        }


        auto ref opIndex(size_t i, size_t j) inout
        in{
            assert(i < rows || rows == 0);
            assert(j < cols || cols == 0);

          static if(hasLength!R)
            assert((mjr == Major.row ? i : j) < _range.length);

          static if(hasLength!(Unqual!(std.range.ElementType!R)))
            assert((mjr == Major.row ? j : i) < _range[i].length);
        }
        body{
          static if(mjr == Major.row)
            return _range[i][j];
          else
            return _range[j][i];
        }

        mixin(defaultExprOps!(false));

      private:
        R _range;
    }
  
    return ToMatrix!()(range);
}
unittest{
    scope(failure) {writefln("Unittest failure :%s(%s)", __FILE__, __LINE__); stdout.flush();}
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}


    auto arr = [[0, 1], [2, 3], [4, 5]];
    auto r1 = toMatrix!(3, 2, Major.row)(arr);
    static assert(isMatrix!(typeof(r1)));
    assert(equal!"equal(a, b)"(r1.toRange, arr));

    auto r2 = arr.toMatrix!(1, 1, Major.row);
    assert(r2[0] == 0);

    auto r3 = arr.toMatrix!(0, 2, Major.row);
    assert(r3.congeal!(3, 2) == r1);

    auto r4 = arr.toMatrix!(2, 0, Major.row);
    assert(equal(r4.congeal!(2, 2)().toFlatten(), [0, 1, 2, 3]));

    auto r5 = arr.toMatrix!(0, 0, Major.row);
    assert(r5 == r1);
}
unittest{
    scope(failure) {writefln("Unittest failure :%s(%s)", __FILE__, __LINE__); stdout.flush();}
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}


    auto arr = [[0, 1], [2, 3], [4, 5]];
    auto r1 = arr.toMatrix!(2, 3, Major.column);
    assert(equal!"equal(a, b)"(r1.transpose.toRange, arr));
    assert(r1[0, 0] == 0); assert(r1[0, 1] == 2); assert(r1[0, 2] == 4);
    assert(r1[1, 0] == 1); assert(r1[1, 1] == 3); assert(r1[1, 2] == 5);

    auto r2 = arr.toMatrix!(1, 1, Major.column);
    assert(r2[0] == 0);

    auto r3 = arr.toMatrix!(wild, 3, Major.column);
    assert(r3 == r1);

    auto r4 = arr.toMatrix!(2, wild, Major.column);
    assert(equal(r4.transpose.toFlatten, [0, 1, 2, 3, 4, 5]));

    auto r5 = arr.toMatrix!(wild, wild, Major.column);
    assert(r5 == r1);
}


/**
単位行列
*/
@property
auto identity(E)()if(isScalar!E)
{
    static struct Identity()
    {
        enum rows = wild;
        enum cols = wild;


        static InferredResult inferSize(size_t i, size_t j)
        {
            if(i == wild && j == wild)
                return InferredResult(false);
            else if(i == wild || j == wild)
                return InferredResult(true, max(i, j), max(i, j));
            else if(i == j)
                return InferredResult(true, i, j);
            else
                return InferredResult(false);
        }


        E opIndex(size_t i, size_t j) inout
        {
            return (i == j) ? (cast(E)1) : (cast(E)0);
        }

        mixin(defaultExprOps!(true));
    }

    return Identity!()();
}
unittest{
    scope(failure) {writefln("Unittest failure :%s(%s)", __FILE__, __LINE__); stdout.flush();}
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}


    auto id = identity!int;

    static assert(typeof(id).inferSize(4, 4).isValid);
    static assert(!typeof(id).inferSize(1, 3).isValid);

    auto m1 = Matrix!(int, 2, 2).init;
    m1.array[] = [0, 1, 2, 3];
    assert(equal((m1 * id).toFlatten, [0, 1, 2, 3]));

    auto id2 = id + id;
    static assert(isMatrix!(typeof(id2)));
    static assert(typeof(id2).inferSize(4, 4).isValid);

    auto id3 = id.congeal!(wild, 2) * id;
    static assert(id3.rows == 2);
    static assert(id3.cols == 2);
    assert(equal(id3.toFlatten, [1, 0, 0, 1]));

    auto ins = id2.congeal!(2, 2);
    static assert(isMatrix!(typeof(ins)));
    assert(equal(ins.toFlatten, [2, 0, 0, 2]));
}


/**
全要素が1な行列を返します。
*/
@property
auto ones(E)()if(isScalar!E)
{
    static struct Ones()
    {
        enum rows = wild;
        enum cols = wild;


        E opIndex(size_t i, size_t j) inout
        {
            return cast(E)1;
        }


        static InferredResult inferSize(size_t i, size_t j)
        {
            if(i == wild && j == wild)
                return InferredResult(false);
            else if(i == wild || j == wild)
                return InferredResult(true, max(i, j), max(i, j));
            else
                return InferredResult(true, i, j);
        }


        //mixin(ExpressionOperators!(ETOSpec.all & ~ETOSpec.matrixMulScalar & ~ETOSpec.scalarMulMatrix).stringMixin);
        mixin(ExpressionOperatorsInferable!(ETOSpec.all & ~ETOSpec.opEquals & ~ETOSpec.toString & ~ETOSpec.matrixMulScalar & ~ETOSpec.scalarMulMatrix).stringMixin);
        const{mixin(ExpressionOperatorsInferable!(ETOSpec.all & ~ETOSpec.matrixMulScalar & ~ETOSpec.scalarMulMatrix).stringMixin);}
        immutable{mixin(ExpressionOperatorsInferable!(ETOSpec.all & ~ETOSpec.opEquals & ~ETOSpec.toString & ~ETOSpec.matrixMulScalar & ~ETOSpec.scalarMulMatrix).stringMixin);}


        auto opBinary(string op : "*", S)(S s) const
        if(isScalar!S)
        {
            return Ns!S(s);
        }


        auto opBinaryRight(string op : "*", S)(S s) const
        if(isScalar!S)
        {
            return Ns!S(s);
        }
    }

    static struct Ns(E)
    {
        enum rows = 0;
        enum cols = 0;


        E opIndex(size_t i, size_t j) inout
        {
            return e;
        }


        static InferredResult inferSize(size_t i, size_t j)
        {
            if(i == wild && j == wild)
                return InferredResult(false);
            else if(i == wild || j == wild)
                return InferredResult(true, max(i, j), max(i, j));
            else
                return InferredResult(true, i, j);
        }


        mixin(ExpressionOperatorsInferable!(ETOSpec.all & ~ETOSpec.opEquals & ~ETOSpec.toString & ~ETOSpec.matrixMulScalar & ~ETOSpec.scalarMulMatrix).stringMixin);
        const{mixin(ExpressionOperatorsInferable!(ETOSpec.all & ~ETOSpec.matrixMulScalar & ~ETOSpec.scalarMulMatrix).stringMixin);}
        immutable{mixin(ExpressionOperatorsInferable!(ETOSpec.all & ~ETOSpec.opEquals & ~ETOSpec.toString & ~ETOSpec.matrixMulScalar & ~ETOSpec.scalarMulMatrix).stringMixin);}


        auto opBinary(string op : "*", S)(S s) const
        if(isScalar!S && is(typeof(s * e)))
        {
            return Ns!(typeof(e * s))(e * s);
        }


        auto opBinaryRight(string op : "*", S)(S s) const
        if(isScalar!S && is(typeof(s * e)))
        {
            return Ns!(typeof(e * s))(s * e);
        }


      private:
        E e;
    }

    static assert(isMatrix!(Ones!()));

    return Ones!().init;
}
unittest{
    scope(failure) {writefln("Unittest failure :%s(%s)", __FILE__, __LINE__); stdout.flush();}
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}


    auto m1 = ones!float;
    assert(m1[0, 1] == 1);

    auto m3 = m1 * 3;
    assert(m3[0, 1] == 3);

    auto m9 = m3 * 3;
    assert(m9[0, 1] == 9);
}


/**
部分行列を返します
*/
auto sub(alias rArray, alias cArray, A)(A mat)
if(isArray!(typeof(rArray)) && isArray!(typeof(cArray)) && isMatrix!A && !isInferableMatrix!A
    && (hasDynamicRows!A || is(typeof({static assert(rArray.find!"a>=b"(A.rows).empty);})))
    && (hasDynamicColumns!A || is(typeof({static assert(cArray.find!"a>=b"(A.cols).empty);}))))
in{
    static if(hasDynamicRows!A)
        assert(rArray.find!"a>=b"(mat.rows).empty);
    static if(hasDynamicColumns!A)
        assert(cArray.find!"a>=b"(mat.cols).empty);
}
body{
  static if(rArray.length == 0 && cArray.length == 0)
    return mat;
  else
  {
    static struct Sub()
    {
      static if(rArray.length == 0)
        alias rows = _mat.rows;
      else
        enum rows = rArray.length;

      static if(cArray.length == 0)
        alias cols = _mat.cols;
      else
        enum cols = cArray.length;


        auto ref opIndex(size_t i, size_t j) inout
        in{
            assert(i < rows);
            assert(j < cols);
        }
        body{
          static if(rArray.length && cArray.length)
            return _mat[rArray[i], cArray[j]];
          else static if(rArray.length)
            return _mat[rArray[i], j];
          else static if(cArray.length)
            return _mat[i, cArray[j]];
          else
            static assert(0);
        }


        mixin(defaultExprOps!(false));


      private:
        A _mat;
    }

    return Sub!()(mat);
  }
}
unittest{
    scope(failure) {writefln("Unittest failure :%s(%s)", __FILE__, __LINE__); stdout.flush();}
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}


    auto m1 = [[0, 1], [2, 3], [4, 5]].toMatrix!(3, 2, Major.row);
    auto s1 = m1.sub!([1, 2], [0]);
    static assert(s1.rows == 2);
    static assert(s1.cols == 1);

    assert(s1[0, 0] == 2);
    assert(s1[1, 0] == 4);


    auto m2 = [[0, 1], [2, 3], [4, 5]].toMatrix!(0, 0, Major.row)();
    auto s2 = sub!((size_t[]).init, (size_t[]).init)(m2);
    assert(m1 == s2.congeal!(3, 2));


    auto m3 = identity!int.congeal!(2, 2)().sub!([0, 0, 0], [0, 0, 0])();
    assert(m3 == ones!int);
}


/**
swizzle : glm参照
*/
auto swizzle(string r, string c, A)(A mat)
if(isMatrix!A)
{
    size_t[] arrGen(string s, size_t size)
    {
        if(s.empty){
            if(size == 0)
                return typeof(return).init;
            else{
                auto arr = new size_t[size];
                foreach(i, ref e; arr)
                    e = i;
                return arr;
            }
        }else{
            auto arr = new size_t[s.length];
            foreach(i, ref e; arr)
                e = s[i] - 'a';
            return arr;
        }
    }

    return sub!(arrGen(r, A.rows), arrGen(c, A.cols))(mat);
}
unittest{
    scope(failure) {writefln("Unittest failure :%s(%s)", __FILE__, __LINE__); stdout.flush();}
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}


    auto a = identity!int.congeal!(2, 2).swizzle!("abab", "baba")();
    assert(equal!"equal(a, b)"(a.toRange, [[0, 1, 0, 1],
                                           [1, 0, 1, 0],
                                           [0, 1, 0, 1],
                                           [1, 0, 1, 0]]));
}


/**
余因子行列
*/


/**
行列の跡(trace)を返します。正方行列についてのみ定義されます
*/
ElementType!A trace(A)(A mat)
if(isMatrix!A && !isInferableMatrix!A && (!(hasStaticRows!A && hasStaticColumns!A) || is(typeof({static assert(A.rows == A.cols);}))))
{
    alias ElementType!A T;
    T sum = cast(T)0;

    foreach(i; 0 .. mat.rows)
        sum += mat[i, i];

    return sum;
}
unittest{
    scope(failure) {writefln("Unittest failure :%s(%s)", __FILE__, __LINE__); stdout.flush();}
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}


    auto m = Matrix!(int, 2, 2)();
    m[0, 0] = 0; m[0, 1] = 1;
    m[1, 0] = 2; m[1, 1] = 3;

    auto tr = m.trace;
    assert(tr == 3);
}


auto dot(V1, V2)(V1 vec1, V2 vec2)
if(isVector!V1 && isVector!V2 && (!(hasStaticLength!V1 && hasStaticLength!V2) || is(typeof({static assert(V1.length == V2.length);}))))
in{
    static if(!(hasStaticLength!V1 && hasStaticLength!V2))
    {
        assert(vec1.length == vec2.length);
    }
}
body{
    alias ElementType!V1 T;
    T sum = cast(T)0;

    foreach(i; 0 .. vec1.length)
        sum += vec1[i] * vec2[i];

    return sum;
}
unittest{
    scope(failure) {writefln("Unittest failure :%s(%s)", __FILE__, __LINE__); stdout.flush();}
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}


    auto rv = Rvector!(int, 3)(),
         cv = Cvector!(int, 3)();

    rv.array[] = [0, 1, 2];
    cv.array[] = [1, 2, 3];

    assert((rv * cv)[0] == 8);  //8
    assert(rv.dot(cv) == 8);
    assert(cv.dot(rv) == 8);

    assert(rv.dot(rv) == 5);
    assert(cv.dot(cv) == 14);
}


/*
auto cross(V1, V2)(V1 vec1, V2 vec2)
if(isVector!V1 && isVector!V2 && )
{
    auto m = vec1.crossAsMatrix(vec2);
}
*/

/**
直積
*/
auto cartesian(V1, V2)(V1 vec1, V2 vec2)
if(isVector!V1 && isVector!V2)
{
    static struct Cartesian()
    {
        alias rows = _vec1.length;
        alias cols = _vec2.length;


        auto opIndex(size_t i, size_t j) const
        in{
            assert(i < rows);
            assert(j < cols);
        }
        body{
            return _vec1[i] * _vec2[j];
        }


        mixin(defaultExprOps!(false));


      private:
        V1 _vec1;
        V2 _vec2;
    }


    return Cartesian!()(vec1, vec2);
}
unittest{
    scope(failure) {writefln("Unittest failure :%s(%s)", __FILE__, __LINE__); stdout.flush();}
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}


    auto v1 = [0, 1, 2, 3].toMatrix!(3, 1);
    auto v2 = [2, 3, 4, 5].toMatrix!(1, 2);

    assert(v1.cartesian(v2) == v1 * v2);
    static assert(hasStaticRows!(typeof(v1.cartesian(v2))));
    static assert(hasStaticColumns!(typeof(v1.cartesian(v2))));
}