/**
Boost::uBlas目指そうと思ったけどuBlasそんなに使ったことないし、
じゃあ独自的な方向でいこうって感じのExpression Templateを使った行列ライブラリ。

特徴は
・ET使ってるから中間オブジェクトの生成がほとんどないらしい
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
opAssign(代入時に演算) : 高速化との兼ね合い。 A = B;を A.opAssign(B);とすると、Bの処理が隠蔽されてしまうので最適化できないのでは？
    A.opAssign -> B.opAssignRight(A)とかA.opOpAssign(B) -> B.opOpAssign(A)とか
内積,外積, 直積,
LU分解, 逆行列,
slice!([r1, r2], [c1, c2]), sub!([r1, r2, r3], [r1, r2, r3, r4]), swizzle
共役転置

帯行列, 対称行列, エルミート行列,
疎行列,
動的大きさを持つ行列(rlength == 0 || clength == 0)
ETを使わない行列

高速化とか
・Blasいくか
    asmかcore.simd.Vectorか外部ライブラリ(C)か。
        外部ライブラリでよくね？自分がやっても遅いし、ユーザーが選択できるし。
            それdbletじゃん…
            じゃあlapackも選択制でいいのでは？version=UserLapackで、lapackを積極的に使うとか
        blasを積極的に使う式テンプレートを使った行列型とか
全要素一括評価：ETの中間オブジェクト云々の良い特徴が消える…？
    A *= 2;だとそんなことは無いけど、ETじゃない。
    A *= B;だと静的大きさでは無理。
        というか外部のBLASを使う際の静的大きさであるメリットはない

Bug:
@safeと@systemの推論のバグ
toString：コンパイルが終わらずメモリ馬鹿食いするバグ

特定行列の最適化(identity * A == AとかA.transpose.transpose == A)


Author: Kazuki Komatsu

License: NYSL
*/

module extml.core;

import std.algorithm,
       std.functional,
       std.range,
       std.traits,
       std.format;

version(unittest) import std.stdio;

/**
デフォルトでの行列が列優先か行優先のどちらでメモリ上に配置されるかをDlinearColumnMajorという
*/
version(ExtmlColumnMajor)
    enum Major defaultMajor = Major.column;
else
    enum Major defaultMajor = Major.row;

//version(unittest) void main(){}


/**
四則演算が定義されている型。
inout(ubyte)からinout(creal)などまで。
*/
template isScalar(T){
    enum bool isScalar = is(typeof(
        {
            T _;
            T a = cast(T)1;
            //a = cast(T)0;
            //a = a;
            
            //a = -a;
            //a = +a;
            
            //a += a;
            //a -= a;
            //a *= a;
            //a /= a;

            //a += 1;
            //a -= 1;
            //a *= 1;
            //a /= 1;

            {auto _t1 = a + a;}
            {auto _t1 = a - a;}
            {auto _t1 = a * a;}
            {auto _t1 = a / a;}

            {auto _t1 = a + 1;}
            {auto _t1 = a - 1;}
            {auto _t1 = a * 1;}
            {auto _t1 = a / 1;}

            //bool b = a < 0;
            //bool b = a == 0;
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


/**
行列型か判定します。
行列型Tは次のコードがコンパイル可能です。

---
T m;

enum rsize = T.rlength;
enum csize = T.clength;

static assert(rsize >= 0);
static assert(csize >= 0);

auto e = m[0, 0];
static assert(isScalar!(typeof(e)));


static if(rsize == 0 || csize == 0)
{
    alias Instance = T.Instantiate!(1, 1);
    enum bool able = Instance.isValid;

    static if(able)
    {
        enum size_t rlength = Instance.rlength;
        enum size_t clength = Instance.clength;
    }
}
---
*/
template isMatrix(T)
{
    enum bool isMatrix = is(typeof({
            T m;

            enum rsize = T.rlength;
            enum csize = T.clength;

            static assert(rsize >= 0);
            static assert(csize >= 0);

            auto e = m[0, 0];
            static assert(isScalar!(typeof(e)));


            static if(rsize == 0 || csize == 0)
            {
                alias Instance = T.Instantiate!(1, 1);
                enum bool able = Instance.isValid;

                static if(able)
                {
                    enum size_t rlength = Instance.rlength;
                    enum size_t clength = Instance.clength;
                }
            }
        }));
}
unittest{
    static struct S
    {
        enum rlength = 1;
        enum clength = 1;

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
        enum rlength = 0;
        enum clength = 0;

        auto opIndex(size_t i, size_t j)
        {
            return i + j;
        }
    }

    static assert(!isMatrix!S);
}
unittest{
    import std.bigint, std.typetuple;
    alias TT = TypeTuple!(ubyte, ushort, uint, ulong,
                           byte,  short,  int,  long,
                          float, double, real,
                          creal, cfloat, cdouble/*,*/
                          /*BigInt*/);

    foreach(T; TT)
    {
        static assert(isMatrix!(Matrix!(T, 3, 3)));
        static assert(isMatrix!(Matrix!(const(T), 3, 3)));
        static assert(isMatrix!(Matrix!(immutable(T), 3, 3)));
    }
}


/**
静的な大きさが特に決まっていない行列をFlexibleMatrixとソースコード上では示しています。
*/
enum FlexibleMatrix{
    row,
    column,
    rowAndColumn,
    any,
}


/**
フレキシブル行列型かどうか判定します。
フレキシブル行列を導入することにより、たとえば単位行列などを抽象的に考えることができ、
さらに、動的な大きさをもつ行列に対して、「どの大きさにもなり得る」という静的な大きさをもたせることができます。

単位行列は、そのもの自身を考える際には、ただ「正方行列」であればよいので、
フレキシブル行列型で表すことができ、次のような実装になります。

Example:
---
struct Identity(T)
{
    enum size_t rlength = 0;    ///行フレキシブル行列
    enum size_t clength = 0;    ///列フレキシブル行列


    ///行列型のプリミティブ
    T opIndex(size_t i, size_t j) const {
        return i == j ? 1 : 0;
    }


    ///フレキシブル行列のプリミティブ
    ///意味は、「単位行列は正方行列でなければいけない」
    ///0の場合はワイルドカードで、推論可能であれば推論するが、推論不可能なら推論しないで0にする。
    static struct Instantiate(size_t i, size_t j)
    {
        static if(i == 0 && j != 0) //iのみがワイルドカード
        {
            enum isValid = true;    //iはjから推論可能
            enum rlength = j;       //rlength = clength = j
            enum clength = j;
        }
        else static if(j == 0 && i != 0) //上と逆
        {
            enum isValid = true;
            enum rlength = i;
            enum clength = i;
        }
        else                        //どちらもワイルドカードでない
        {
            enum isValid = i == j;  //i == jなら正方行列なので正当

            static if(isValid)
            {
                enum rlength = i;   //正当な場合にはrlengthとclengthを宣言
                enum clength = j;
            }
        }
    }
}
---


列の大きさが動的な大きさを持つ行列は次のように表せます。

Example:
---
struct DynamicMatrix(T, size_t rsize)
{
    enum size_t rlength = rsize;
    enum size_t clength = 0;


    T opIndex(size_t i, size_t j) inout 
    in{
        assert(i < 3);
        assert(3 * j + i < _mem.length);
    }
    body{
        return _mem[3 * j + i];
    }


    static struct Instantiate(size_t r, size_t c)
    {
        static if(r == 0 || r == rsize){
            enum isValid = true;
            enum rlength = r;
            enum clength = c;
        }else
            enum isValid = false;
    }

  private:
    T[]
}
---
*/
template isFlexibleMatrix(T, FlexibleMatrix Type = FlexibleMatrix.any)
{
    static if(Type == FlexibleMatrix.any){
        enum bool isFlexibleMatrix = isMatrix!T
                                  && is(typeof({
                                        static assert(T.rlength == 0 || T.clength == 0);
                                        alias Instance = T.Instantiate!(1, 1);
                                        enum bool able = Instance.isValid;

                                        static if(able)
                                        {
                                            enum size_t rlength = Instance.rlength;
                                            enum size_t clength = Instance.clength;
                                        }
                                    }));
    }
    else
    {
        static if(Type == FlexibleMatrix.rowAndColumn)
            enum bool _isFlexibleMatrix = isMatrix!T && T.rlength == 0 && T.clength == 0;
        else static if(Type == FlexibleMatrix.row)
            enum bool _isFlexibleMatrix = isMatrix!T && T.rlength == 0;
        else static if(Type == FlexibleMatrix.column)
            enum bool _isFlexibleMatrix = isMatrix!T && T.clength == 0;

        enum bool isFlexibleMatrix = _isFlexibleMatrix && is(typeof({
                alias Instance = T.Instantiate!(1, 1);
                enum bool able = Instance.isValid;
                static if(able){
                    enum size_t rlength = Instance.rlength;
                    enum size_t clength = Instance.clength;
                }
            }));
    }
}
unittest{
    static struct S(size_t r, size_t c)
    {
        enum rlength = r;
        enum clength = c;

        int opIndex(size_t, size_t){
            return 1;
        }

        static struct Instantiate(size_t i, size_t j)
        {
            static if(i == 0 && j != 0)
            {
                enum isValid = true;
                enum rlength = j;
                enum clength = j;
            }
            else static if(j == 0 && i != 0)
            {
                enum isValid = true;
                enum rlength = i;
                enum clength = i;
            }
            else
            {
                enum isValid = i == j;
                
                static if(isValid)
                {
                    enum rlength = i;
                    enum clength = j;
                }
            }
        }
    }

    static assert(isFlexibleMatrix!(S!(0, 0), FlexibleMatrix.rowAndColumn));
    static assert(isFlexibleMatrix!(S!(0, 1), FlexibleMatrix.row));
    static assert(isFlexibleMatrix!(S!(1, 0), FlexibleMatrix.column));

    static assert(!isFlexibleMatrix!(S!(1, 1), FlexibleMatrix.row));
    static assert(!isFlexibleMatrix!(S!(1, 1), FlexibleMatrix.column));
    static assert(!isFlexibleMatrix!(S!(1, 1), FlexibleMatrix.rowAndColumn));
}
unittest{
    static assert(isFlexibleMatrix!(typeof(identity!int)));
}


/**
フレキシブル行列のある大きさに対するインスタンスを生成します
*/
auto instantiate(size_t r, size_t c, A)(A fmat)
if(isFlexibleMatrix!(A, FlexibleMatrix.any) && A.Instantiate!(r, c).isValid)
{
    ///ditto
    static struct FlexibleMatrixInstance()
    {
        A _fmat;

        enum rlength = r;
        enum clength = c;

        auto ref opIndex(size_t i, size_t j)
        in{
            assert(i < rlength || rlength == 0);
            assert(j < clength || clength == 0);
        }
        body{
            return _fmat[i, j];
        }


        auto ref opIndex(size_t i, size_t j) const 
        in{
            assert(i < rlength || rlength == 0);
            assert(j < clength || clength == 0);
        }
        body{
            return _fmat[i, j];
        }


        auto ref opIndex(size_t i, size_t j) immutable
        in{
            assert(i < rlength || rlength == 0);
            assert(j < clength || clength == 0);
        }
        body{
            return _fmat[i, j];
        }


      mixin ExpressionTemplateOperators!(ETOSpec.all);

      static if(is(typeof(_opIndex)))
        mixin(_opIndex);


      static if(this.rlength == 0 || this.clength == 0)
      {
        static struct Instantiate(size_t R, size_t C)
        {
          private:
            enum bool _isValid = isEqOrEitherEq0(FlexibleMatrixInstance.rlength, R)
                              && isEqOrEitherEq0(FlexibleMatrixInstance.clength, C);

            static if(_isValid)
            {
                enum size_t _rlength = max(FlexibleMatrixInstance.rlength, R);
                enum size_t _clength = max(FlexibleMatrixInstance.clength, C);
            }

          public:
            static if(_isValid)
                enum bool isValid = A.Instantiate!(_rlength, _clength).isValid;
            else
                enum bool isValid = false;

            static if(isValid)
            {
                enum size_t rlength = A.Instantiate!(_rlength, _clength).rlength;
                enum size_t clength = A.Instantiate!(_rlength, _clength).clength;
            }

        }
      }
    }

    return FlexibleMatrixInstance!()(fmat);
}

unittest{
    auto id = identity!int.instantiate!(3, 3);
    static assert(id.rlength == 3);
    static assert(id.clength == 3);
    static assert(isMatrix!(typeof(id)));

    assert(equal!"equal(a, b)"(id.toRange, [[1, 0, 0], [0, 1, 0], [0, 0, 1]]));
}


/**
ベクトル型かどうか判定します。
*/
template isVector(V)
{
    enum isVector = isMatrix!V
                 && is(typeof({
                        static assert(V.rlength == 1 || V.clength == 1);

                        enum len = V.length;
                        static assert(len == V.rlength * V.clength);

                        V v;

                        auto a = v[0, 0];
                        auto b = v[0];
                        static assert(is(typeof(a) == typeof(b)));
                    }));
}
unittest{
    static assert(isVector!(Rvector1ui));
    static assert(isVector!(Rvector3ui));
    static assert(isVector!(Cvector1ui));
    static assert(isVector!(Cvector3ui));
}


/**
行列型の要素の型を取得します。

Example:
---
static assert(is(ElementType!(NMatrix2S!int) == int));
---
*/
template ElementType(A)if(isMatrix!A)
{
    alias typeof(A.init[0, 0]) ElementType;
}
unittest{
    static struct S
    {
        enum rlength = 1;
        enum clength = 1;

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
    static assert(hasLvalueElements!(Matrix3i));
}
unittest{
    static struct M
    {
        enum rlength = 1;
        enum clength = 1;

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
    static assert(hasAssignableElements!(Matrix3i));
}
unittest{
    static struct M
    {
        enum rlength = 1;
        enum clength = 1;

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
その2項演算子が正しいかどうかを判定します。

Example:
---
unittest{
    static assert(isValidOperator!(Matrix1x1i, "*", Matrix1x1f));
    static assert(isValidOperator!(Matrix2x1i, "*", Matrix1x2f));
    static assert(isValidOperator!(Matrix1x2i, "*", Matrix2x1f));
    static assert(!isValidOperator!(Matrix1x2i, "*", Matrix1x2f));
}
unittest{
    static assert(isValidOperator!(Matrix1x1i, "+", Matrix1x1f));
    static assert(!isValidOperator!(Matrix2x1i, "+", Matrix1x2f));
    static assert(!isValidOperator!(Matrix1x2i, "+", Matrix2x1f));
    static assert(isValidOperator!(Matrix1x2i, "+", Matrix1x2f));
}
unittest{
    static assert(isValidOperator!(Matrix1x1i, "-", Matrix1x1f));
    static assert(!isValidOperator!(Matrix2x1i, "-", Matrix1x2f));
    static assert(!isValidOperator!(Matrix1x2i, "-", Matrix2x1f));
    static assert(isValidOperator!(Matrix1x2i, "-", Matrix1x2f));
}
unittest{
    static assert(!isValidOperator!(int, "-", Matrix1x1f));
    static assert(!isValidOperator!(int, "+", Matrix1x2f));
    static assert(isValidOperator!(int, "*", Matrix2x1f));
    static assert(!isValidOperator!(Matrix2x1i, "-", creal));
    static assert(!isValidOperator!(Matrix1x2i, "+", creal));
    static assert(isValidOperator!(Matrix1x2i, "*", creal));
}
---
*/
template isValidOperator(L, string op, R)
{
    static assert(isScalar!L || isMatrix!L);
    static assert(isScalar!R || isMatrix!R);
    static assert(isMatrix!L || isMatrix!R);

    static if(isMatrix!L)
    {
        static if(isMatrix!R)
        {
            static if(isFlexibleMatrix!(L) && isFlexibleMatrix!(R))
            {
                static if(op == "+" || op == "-")
                    enum bool _isValidOperator = is(typeof({
                                                        alias LI1 = L.Instantiate!(R.rlength, R.clength);
                                                        static assert(LI1.isValid);

                                                        alias RI1 = R.Instantiate!(LI1.rlength, LI1.clength);
                                                        static assert(RI1.isValid);

                                                        alias LI2 = L.Instantiate!(RI1.rlength, RI1.clength);
                                                        static assert(LI2.isValid);
                                                    }));
                else static if(op == "*")
                {
                    static if(L.clength == 0 && R.rlength == 0)
                    {
                        //pragma(msg, "ColumnFlexibleMatrix * RowFlexibleMatrix() cannot instantiate.");
                        //pragma(msg, "Please use FlexibleMatrix.instantiate!(r, c)");
                        enum bool _isValidOperator = false;
                    }
                    else static if(L.clength == 0)
                    {
                        enum bool _isValidOperator = L.Instantiate!(0, R.rlength).isValid;
                    }
                    else static if(R.rlength == 0)
                    {
                        enum bool _isValidOperator = R.Instantiate!(L.clength, 0).isValid;
                    }
                    else
                    {
                        //pragma(msg, "Invalid case.");
                        static assert(0);
                    }
                }
                else
                    enum bool _isValidOperator = false;
            }
            else static if(isFlexibleMatrix!(L))
            {
                static if(op == "+" || op == "-")
                    enum bool _isValidOperator = L.Instantiate!(R.rlength, R.clength).isValid
                                              && is(typeof({
                                                    static assert(isEqOrEitherEq0(L.Instantiate!(R.rlength, R.clength).rlength, R.rlength));
                                                    static assert(isEqOrEitherEq0(L.Instantiate!(R.rlength, R.clength).clength, R.clength));
                                                }));
                else static if(op == "*")
                    enum bool _isValidOperator = L.Instantiate!(0, R.rlength).isValid
                                             && is(typeof({
                                                    static assert(isEqOrEitherEq0(L.Instantiate!(0, R.rlength).clength, R.rlength));
                                                }));
                else
                    enum bool _isValidOperator = false;
            }
            else static if(isFlexibleMatrix!(R))
            {
                static if(op == "+" || op == "-")
                    enum bool _isValidOperator = R.Instantiate!(L.rlength, L.clength).isValid
                                             && is(typeof({
                                                    static assert(isEqOrEitherEq0(R.Instantiate!(L.rlength, L.clength).rlength, L.rlength));
                                                    static assert(isEqOrEitherEq0(R.Instantiate!(L.rlength, L.clength).clength, L.clength));
                                                }));
                else static if(op == "*")
                    enum bool _isValidOperator = R.Instantiate!(L.clength, 0).isValid
                                              && is(typeof({
                                                    static assert(isEqOrEitherEq0(R.Instantiate!(L.clength, 0).rlength, L.clength));
                                                }));
                else
                    enum bool _isValidOperator = false;
            }
            else
            {
                static if(op == "+" || op == "-")
                    enum bool _isValidOperator = L.rlength == R.rlength
                                              && L.clength == R.clength;
                else static if(op == "*")
                    enum bool _isValidOperator = L.clength == R.rlength;
                else
                    enum bool _isValidOperator = false;
            }

            enum bool isValidOperator = is(typeof(mixin("L.init[0, 0] " ~ op ~ "R.init[0, 0]")))
                                     && _isValidOperator;
        }else
            enum bool isValidOperator = op == "*" && is(typeof(mixin("L.init[0, 0] " ~ op ~ "R.init")));
    }
    else
    { // isScalar
        //static assert(isMatrix!R);
        enum bool isValidOperator = op == "*" && is(typeof(mixin("L.init " ~ op ~ " R.init[0, 0]")));
    }
}
unittest{
    static assert(isValidOperator!(Matrix1x1i,  "*", Matrix1x1f));
    static assert(isValidOperator!(Matrix2x1i,  "*", Matrix1x2f));
    static assert(isValidOperator!(Matrix1x2i,  "*", Matrix2x1f));
    static assert(!isValidOperator!(Matrix1x2i, "*", Matrix1x2f));
}
unittest{
    static assert(isValidOperator!(Matrix1x1i,  "+", Matrix1x1f));
    static assert(!isValidOperator!(Matrix2x1i, "+", Matrix1x2f));
    static assert(!isValidOperator!(Matrix1x2i, "+", Matrix2x1f));
    static assert(isValidOperator!(Matrix1x2i,  "+", Matrix1x2f));
}
unittest{
    static assert(isValidOperator!(Matrix1x1i,  "-", Matrix1x1f));
    static assert(!isValidOperator!(Matrix2x1i, "-", Matrix1x2f));
    static assert(!isValidOperator!(Matrix1x2i, "-", Matrix2x1f));
    static assert(isValidOperator!(Matrix1x2i,  "-", Matrix1x2f));
}
unittest{
    static assert(!isValidOperator!(int, "-", Matrix1x1f));
    static assert(!isValidOperator!(int, "+", Matrix1x2f));
    static assert(isValidOperator!(int,  "*", Matrix2x1f));
    static assert(!isValidOperator!(Matrix2x1i, "-", creal));
    static assert(!isValidOperator!(Matrix1x2i, "+", creal));
    static assert(isValidOperator!(Matrix1x2i,  "*", creal));
}
unittest{
    static assert(isValidOperator!(typeof(identity!int),  "+", Matrix1x1f));
    static assert(!isValidOperator!(typeof(identity!int), "+", Matrix1x2f));
    static assert(!isValidOperator!(typeof(identity!int), "+", Matrix2x1f));
    static assert(isValidOperator!(typeof(identity!int),  "+", Matrix2x2f));

    static assert(isValidOperator!(typeof(identity!int),  "-", Matrix1x1f));
    static assert(!isValidOperator!(typeof(identity!int), "-", Matrix1x2f));
    static assert(!isValidOperator!(typeof(identity!int), "-", Matrix2x1f));
    static assert(isValidOperator!(typeof(identity!int),  "-", Matrix2x2f));

    static assert(isValidOperator!(typeof(identity!int),  "*", Matrix1x1f));
    static assert(isValidOperator!(typeof(identity!int),  "*", Matrix1x2f));
    static assert(isValidOperator!(typeof(identity!int),  "*", Matrix2x1f));
    static assert(isValidOperator!(typeof(identity!int),  "*", Matrix2x2f));
}


/**
二項演算子での結果を調べます。

Example:
---
alias Result!(Matrix2x1i, "*", Matrix1x2f) R;
static assert(R.rlength == 1);
static assert(R.clength == 1);
---
*/
template ResultSize(L, string op, R)if(isValidOperator!(L, op, R))
{
    static assert(isScalar!L || isMatrix!L);
    static assert(isScalar!R || isMatrix!R);
    static assert(isMatrix!L || isMatrix!R);

    static if(isMatrix!L)
    {
        static if(isMatrix!R)
        {
            static if(isFlexibleMatrix!L && isFlexibleMatrix!R)
            {
                static if(op == "+" || op == "-")
                {
                    enum rlength = max(L.rlength, R.rlength);
                    enum clength = max(L.clength, R.clength);
                }
                else static if(op == "*")
                {
                    /*
                    static if(L.clength == 0 && R.rlength == 0)
                    {}
                    else */static if(L.clength == 0)
                    {
                        enum rlength = L.Instantiate!(0, R.rlength).rlength;
                        enum clength = R.Instantiate!(R.rlength, 0).clength;
                    }
                    else static if(R.rlength == 0)
                    {
                        enum rlength = L.Instantiate!(0, L.clength).rlength;
                        enum clength = R.Instantiate!(L.clength, 0).clength;
                    }
                }
            }
            else static if(isFlexibleMatrix!L)
            {
                static if(op == "+" || op == "-")
                {
                    enum rlength = R.rlength;
                    enum clength = R.clength;
                }
                else static if(op == "*")
                {
                    enum rlength = L.Instantiate!(0, R.rlength).rlength;
                    enum clength = R.clength;
                }
            }
            else static if(isFlexibleMatrix!R)
            {
                static if(op == "+" || op == "-")
                {
                    enum rlength = L.rlength;
                    enum clength = L.clength;
                }
                else static if(op == "*")
                {
                    enum rlength = L.rlength;
                    enum clength = R.Instantiate!(L.clength, 0).clength;
                }
            }
            else
            {
                enum rlength = L.rlength;
                enum clength = R.clength;
            }
        }
        else
        {
            alias typeof(mixin("L.init[0, 0]" ~ op ~ "R.init")) ElementType;
            enum rlength = L.rlength;
            enum clength = L.clength;
        }
    }
    else
    {
        alias typeof(mixin("L.init" ~ op ~ "R.init[0, 0]")) ElementType;
        enum rlength = R.rlength;
        enum clength = R.clength;
    }
}
unittest{
    alias ResultSize!(Matrix1x2i, "*", Matrix2x1i) R;
    static assert(R.rlength == 1);
    static assert(R.clength == 1);
}
unittest{
    alias ResultSize!(Matrix2x1i, "*", Matrix1x2f) R;
    static assert(R.rlength == 2);
    static assert(R.clength == 2);
}
unittest{
    alias I = typeof(identity!int);
    alias ResultSize!(I, "+", I) R;
    static assert(R.rlength == 0);
    static assert(R.clength == 0);
}
unittest{
    alias ResultSize!(typeof(identity!int.instantiate!(0, 2)), "*", typeof(identity!int)) R;
    static assert(R.rlength == 2);
    static assert(R.clength == 2);
}


/**
Expression Template Operator Species : 式テンプレート演算子の種類
*/
enum ETOSpec{
    matrixAddMatrix = (1 << 0),
    matrixSubMatrix = (1 << 1),
    matrixMulMatrix = (1 << 2),
    matrixMulScalar = (1 << 3),
    scalarMulMatrix = (1 << 4),
    matrixEquMatrix = (1 << 5),
    all = (1 << 6) -1,
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
    static if(op == "+")
        enum ETOSpec ETOperatorSpec = ETOSpec.matrixAddMatrix;
    else static if(op == "-")
        enum ETOSpec ETOperatorSpec = ETOSpec.matrixSubMatrix;
    else static if(op == "*")
    {
        static if(isMatrix!A)
        {
            static if(isMatrix!B)
                enum ETOSpec ETOperatorSpec = ETOSpec.matrixMulMatrix;
            else
                enum ETOSpec ETOperatorSpec = ETOSpec.matrixMulScalar;
        }else
            enum ETOSpec ETOperatorSpec = ETOSpec.scalarMulMatrix;
    }
}
unittest{
    static assert(ETOperatorSpec!(Matrix2i, "+", Matrix2i) == ETOSpec.matrixAddMatrix);
    static assert(ETOperatorSpec!(Matrix2i, "-", Matrix2i) == ETOSpec.matrixSubMatrix);
    static assert(ETOperatorSpec!(Matrix2i, "*", Matrix2i) == ETOSpec.matrixMulMatrix);
    static assert(ETOperatorSpec!(Matrix2i, "*", int) == ETOSpec.matrixMulScalar);
    static assert(ETOperatorSpec!(int, "*", Matrix1x2i) == ETOSpec.scalarMulMatrix);
}



///ditto
struct MatrixOperator(string s, A, B)
{
    enum size_t rlength = ResultSize!(A, s, B).rlength;
    enum size_t clength = ResultSize!(A, s, B).clength;
    enum ETOSpec etoSpec = ETOperatorSpec!(A, s, B);

    A lhs;
    B rhs;


    auto ref opIndex(size_t i, size_t j)
    in{
        assert(i < rlength || rlength == 0);
        assert(j < clength || clength == 0);
    }
    body{
      static if(etoSpec == ETOSpec.matrixAddMatrix)
        return this.lhs[i, j] + this.rhs[i, j];
      else static if(etoSpec == ETOSpec.matrixSubMatrix)
        return this.lhs[i, j] - this.rhs[i, j];
      else static if(etoSpec == ETOSpec.matrixMulMatrix)
      {
        Unqual!(typeof(this.lhs[0, 0] * this.rhs[0, 0])) sum = cast(Unqual!(typeof(this.lhs[0, 0] * this.rhs[0, 0])))0;
        foreach(k; 0 .. A.clength)
            sum += this.lhs[i, k] * this.rhs[k, j];
        return sum;
      }else static if(etoSpec == ETOSpec.matrixMulScalar)
        return this.lhs[i, j] * this.rhs;
      else static if(etoSpec == ETOSpec.scalarMulMatrix)
        return this.lhs * this.rhs[i, j];
      else
        static assert(0);
    }


    auto ref opIndex(size_t i, size_t j) const 
    in{
        assert(i < rlength || rlength == 0);
        assert(j < clength || clength == 0);
    }
    body{
      static if(etoSpec == ETOSpec.matrixAddMatrix)
        return this.lhs[i, j] + this.rhs[i, j];
      else static if(etoSpec == ETOSpec.matrixSubMatrix)
        return this.lhs[i, j] - this.rhs[i, j];
      else static if(etoSpec == ETOSpec.matrixMulMatrix)
      {
        Unqual!(typeof(this.lhs[0, 0] * this.rhs[0, 0])) sum = cast(Unqual!(typeof(this.lhs[0, 0] * this.rhs[0, 0])))0;
        foreach(k; 0 .. A.clength)
            sum += this.lhs[i, k] * this.rhs[k, j];
        return sum;
      }else static if(etoSpec == ETOSpec.matrixMulScalar)
        return this.lhs[i, j] * this.rhs;
      else static if(etoSpec == ETOSpec.scalarMulMatrix)
        return this.lhs * this.rhs[i, j];
      else
        static assert(0);
    }


    auto ref opIndex(size_t i, size_t j) immutable
    in{
        assert(i < rlength || rlength == 0);
        assert(j < clength || clength == 0);
    }
    body{
      static if(etoSpec == ETOSpec.matrixAddMatrix)
        return this.lhs[i, j] + this.rhs[i, j];
      else static if(etoSpec == ETOSpec.matrixSubMatrix)
        return this.lhs[i, j] - this.rhs[i, j];
      else static if(etoSpec == ETOSpec.matrixMulMatrix)
      {
        Unqual!(typeof(this.lhs[0, 0] * this.rhs[0, 0])) sum = cast(Unqual!(typeof(this.lhs[0, 0] * this.rhs[0, 0])))0;
        foreach(k; 0 .. A.clength)
            sum += this.lhs[i, k] * this.rhs[k, j];
        return sum;
      }else static if(etoSpec == ETOSpec.matrixMulScalar)
        return this.lhs[i, j] * this.rhs;
      else static if(etoSpec == ETOSpec.scalarMulMatrix)
        return this.lhs * this.rhs[i, j];
      else
        static assert(0);
    }



  static if(this.rlength == 0 || this.clength == 0)
  {
    static struct Instantiate(size_t r, size_t c)
    if(etoSpec == ETOSpec.matrixAddMatrix || etoSpec == ETOSpec.matrixSubMatrix)
    {
        static assert(isFlexibleMatrix!A && isFlexibleMatrix!B);    //どっちもFM
      private:
        alias AI = A.Instantiate!(r, c);
        alias BI = B.Instantiate!(r, c);

      public:
        enum bool isValid = AI.isValid && BI.isValid
                           && isEqOrEitherEq0(AI.rlength, BI.rlength)
                           && isEqOrEitherEq0(AI.clength, BI.clength);

        static if(isValid)
        {
            enum size_t rlength = AI.rlength;
            enum size_t clength = AI.clength;
        }
    }


    static struct Instantiate(size_t r, size_t c)
    if(etoSpec == ETOSpec.matrixMulMatrix)
    {
        static assert(isFlexibleMatrix!A || isFlexibleMatrix!B);    //どちらかがFM

        static if(isFlexibleMatrix!A && isFlexibleMatrix!B)
        {
            private
            enum bool _isValid = A.Instantiate!(r, c).isValid
                             && B.Instantiate!(r, c).isValid;

            static if(_isValid)
                enum bool isValid = _isValid && isEqOrEitherEq0(A.Instantiate!(r, c).clength, B.Instantiate!(r, c).rlength);
            else
                enum bool isValid = false;

            static if(isValid)
            {
                enum rlength = A.Instantiate!(r, c).rlength;
                enum clength = B.Instantiate!(r, c).clength;
            }
        }
        else static if(isFlexibleMatrix!A)
        {
            private
            enum bool _isValid = A.Instantiate!(r, c).isValid;

            static if(_isValid)
                enum bool isValid = _isValid && isEqOrEitherEq0(A.Instantiate!(r, c).clength, B.rlength);
            else
                enum bool isValid = false;

            static if(isValid)
            {
                enum rlength = A.Instantiate!(r, c).rlength;
                enum clength = B.clength;
            }
        }
        else static if(isFlexibleMatrix!B)
        {
            private
            enum bool _isValid = B.Instantiate!(r, c).isValid;

            static if(_isValid)
                enum bool isValid = _isValid && isEqOrEitherEq0(A.clength, B.Instantiate!(r, c).rlength);
            else
                enum bool isValid = false;

            static if(isValid)
            {
                enum rlength = A.rlength;
                enum clength = B.Instantiate!(r, c).clength;
            }
        }
        else
            static assert(0);
    }


    static struct Instantiate(size_t r, size_t c)
    if(etoSpec == ETOSpec.matrixMulScalar)
    {
        static assert(isFlexibleMatrix!A);

        enum isValid = A.Instantiate!(r, c).isValid;
        static if(isValid)
        {
            enum rlength = A.Instantiate!(r, c).rlength;
            enum clength = A.Instantiate!(r, c).clength;
        }
    }


    static struct Instantiate(size_t r, size_t c)
    if(etoSpec == ETOSpec.scalarMulMatrix)
    {
        static assert(isFlexibleMatrix!B);

        enum isValid = B.Instantiate!(r, c).isValid;

        static if(isValid)
        {
            enum rlength = B.Instantiate!(r, c).rlength;
            enum clength = B.Instantiate!(r, c).clength;
        }
    }
  }


    mixin ExpressionTemplateOperators!(ETOSpec.all);

  static if(is(typeof(_opIndex)))
        mixin(_opIndex);
}


/**
行列に関する演算子の式
*/
template matrixOperator(string s, A, B)
if(isValidOperator!(A, s, B))
{
    MatrixOperator!(s, A, B) matrixOperator(A a, B b)
    {
        return MatrixOperator!(s, A, B)(a, b);
    }

    //static assert(isMatrix!(MatrixOperator!(s, A, B)));
}


unittest{
    scope(failure) writefln("Unittest failure : ", __FILE__, " : ", __LINE__);
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}

    auto m1 = [0, 1, 2, 3].toMatrix!(2, 2),
         m2 = [2, 4, 1, 3].toMatrix!(2, 2);

    assert(equal((m1 + m2).toFlatten, [2, 5, 3, 6]));
    assert(equal((m1 - m2).toFlatten, [-2, -3, 1, 0]));
    assert(equal((m1 * m2).toFlatten, [1, 3, 7, 17]));
    assert(equal((m1 *  2).toFlatten, [0, 2, 4, 6]));
    assert(equal(( 3 * m2).toFlatten, [6, 12, 3, 9]));
}
unittest{
    scope(failure) writefln("Unittest failure : ", __FILE__, " : ", __LINE__);
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}


    auto m1 = [0, 1, 2, 3].toMatrix!(2, 2);

    auto m = (m1 + m1) * (m1 + m1 * 3);
    assert(m.toFlatten.equal([16, 24, 48, 88]));
}


/**
specにいれた種類の演算子をオーバーロードします。

Example:
---
struct Identity(T){
    enum rlength = 0;
    enum clength = 0;

    T opIndex(size_t i, size_t j){
        return i == j ? cast(T)0 : cast(T)1;
    }

    mixin ExpressionTemplateOperators!(ETOSpec.all);
}
---
*/
mixin template ExpressionTemplateOperators(ETOSpec spec)
{
  static if(spec & ETOSpec.matrixAddMatrix)
  {
    auto opBinary(string op : "+", Rhs)(Rhs mat)
    if(isMatrix!Rhs)
    {
        static assert(isValidOperator!(typeof(this), op, Rhs));
        return matrixOperator!"+"(this, mat);
    }

    //inoutだとconstやimmutableが伝搬できない
    auto opBinary(string op : "+", Rhs)(const Rhs mat)
    if(isMatrix!Rhs)
    {
        static assert(isValidOperator!(typeof(this), op, const(Rhs)));
        return matrixOperator!"+"(this, mat);
    }


    auto opBinary(string op : "+", Rhs)(immutable Rhs mat)
    if(isMatrix!Rhs)
    {
        static assert(isValidOperator!(typeof(this), op, immutable(Rhs)));
        return matrixOperator!"+"(this, mat);
    }
  }


  static if(spec & ETOSpec.matrixSubMatrix)
  {
    auto opBinary(string op : "-", Rhs)(Rhs mat)
    if(isMatrix!Rhs)
    {
        static assert(isValidOperator!(typeof(this), op, Rhs));
        return matrixOperator!"-"(this, mat);
    }


    auto opBinary(string op : "-", Rhs)(const Rhs mat)
    if(isMatrix!Rhs)
    {
        static assert(isValidOperator!(typeof(this), op, const(Rhs)));
        return matrixOperator!"-"(this, mat);
    }


    auto opBinary(string op : "-", Rhs)(immutable Rhs mat)
    if(isMatrix!Rhs)
    {
        static assert(isValidOperator!(typeof(this), op, immutable(Rhs)));
        return matrixOperator!"-"(this, mat);
    }
  }


  static if(spec & ETOSpec.matrixMulMatrix)
  {
    auto opBinary(string op : "*", Rhs)(Rhs mat)
    if(isMatrix!Rhs)
    {
        static assert(isValidOperator!(typeof(this), op, Rhs));
        return matrixOperator!"*"(this, mat);
    }


    auto opBinary(string op : "*", Rhs)(const Rhs mat)
    if(isMatrix!Rhs)
    {
        static assert(isValidOperator!(typeof(this), op, const(Rhs)));
        return matrixOperator!"*"(this, mat);
    }


    auto opBinary(string op : "*", Rhs)(immutable Rhs mat)
    if(isMatrix!Rhs)
    {
        static assert(isValidOperator!(typeof(this), op, immutable(Rhs)));
        return matrixOperator!"*"(this, mat);
    }
  }


  static if(spec & ETOSpec.matrixMulScalar)
  {
    auto opBinary(string op : "*", S)(S s)
    if(isScalar!S)
    {
        static assert(isValidOperator!(typeof(this), op, S));
        return matrixOperator!"*"(this, s);
    }


    auto opBinary(string op : "*", S)(const S s)
    if(isScalar!S)
    {
        static assert(isValidOperator!(typeof(this), op, const(S)));
        return matrixOperator!"*"(this, s);
    }


    auto opBinary(string op : "*", S)(immutable S s)
    if(isScalar!S)
    {
        static assert(isValidOperator!(typeof(this), op, immutable(S)));
        return matrixOperator!"*"(this, s);
    }
  }


  static if(spec & ETOSpec.scalarMulMatrix)
  {
    auto opBinaryRight(string op : "*", S)(S s)
    if(isScalar!S)
    {
        static assert(isValidOperator!(S, op, typeof(this)));
        return matrixOperator!"*"(s, this);
    }


    auto opBinaryRight(string op : "*", S)(const S s)
    if(isScalar!S)
    {
        static assert(isValidOperator!(typeof(this), op, const(S)));
        return matrixOperator!"*"(s, this);
    }


    auto opBinaryRight(string op : "*", S)(immutable S s)
    if(isScalar!S)
    {
        static assert(isValidOperator!(typeof(this), op, immutable(S)));
        return matrixOperator!"*"(s, this);
    }
  }

  static if(spec & ETOSpec.matrixEquMatrix)
  {
    //pragma(msg, __FILE__, "(", __LINE__, "): @system");
    bool opEquals(Rhs)(ref const Rhs mat) @system
    if(isMatrix!Rhs)
    {
        static assert(isValidOperator!(typeof(this), "+", const(Rhs))
                   && ResultSize!(typeof(this), "+", const(Rhs)).rlength != 0
                   && ResultSize!(typeof(this), "+", const(Rhs)).clength != 0);

        foreach(i; 0 .. Rhs.rlength)
            foreach(j; 0 .. Rhs.clength)
                if(this[i, j] != mat[i, j])
                    return false;
        return true;
    }


    //pragma(msg, __FILE__, "(", __LINE__, "): @system");
    bool opEquals(Rhs)(const Rhs mat) @system
    if(isMatrix!Rhs)
    {
        static assert(isValidOperator!(typeof(this), "+", const(Rhs))
                   && ResultSize!(typeof(this), "+", const(Rhs)).rlength != 0
                   && ResultSize!(typeof(this), "+", const(Rhs)).clength != 0);

        foreach(i; 0 .. max(Rhs.rlength, rlength))
            foreach(j; 0 .. max(Rhs.clength, clength))
                if(this[i, j] != mat[i, j])
                    return false;
        return true;
    }
  }

  static if(typeof(this).rlength == 1 && typeof(this).clength == 1)
  {
    @property auto ref _scalarValue(){return this[0, 0];}

    @property auto ref _scalarValue() const {return this[0, 0];}

    @property auto ref _scalarValue() immutable {return this[0, 0];}

    alias _scalarValue this;

    alias rlength length;
    alias length opDollar;

    enum _opIndex = q{  //dmd bug : オーバーロードできない
        auto ref opIndex(size_t i)
        in{
            assert(i == 0);
        }
        body{
            return this[0, 0];
        }


        auto ref opIndex(size_t i) const
        in{
            assert(i == 0);
        }
        body{
            return this[0, 0];
        }


        auto ref opIndex(size_t i) immutable
        in{
            assert(i == 0);
        }
        body{
            return this[0, 0];
        }
    };
  }
  else static if(typeof(this).rlength == 1)
  {
    enum _opIndex = q{  //dmd bug
        auto ref opIndex(size_t i)
        in{
            assert(i < clength || clength == 0);
        }
        body{
            return this[0, i];
        }


        auto ref opIndex(size_t i) const
        in{
            assert(i < clength || clength == 0);
        }
        body{
            return this[0, i];
        }


        auto ref opIndex(size_t i) immutable
        in{
            assert(i < clength || clength == 0);
        }
        body{
            return this[0, i];
        }
    };


    alias clength length;
    alias length opDollar;
  }
  else static if(typeof(this).clength == 1)
  {
    enum _opIndex = q{  //dmd bug
        auto ref opIndex(size_t i)
        in{
            assert(i < rlength || rlength == 0);
        }
        body{
            return this[i, 0];
        }


        auto ref opIndex(size_t i) const
        in{
            assert(i < rlength || rlength == 0);
        }
        body{
            return this[i, 0];
        }


        auto ref opIndex(size_t i) immutable
        in{
            assert(i < rlength || rlength == 0);
        }
        body{
            return this[i, 0];
        }
    };


    alias rlength length;
    alias length opDollar;
  }
  else
  {
    //pragma(msg, __FILE__, "(", __LINE__, "): Vector!(E) opIndex(size_t i)...");
  }

/*   //dmd bug : toStringをONにすると、メモリをバカ食いする事象*/
    @property
    void toString(scope void delegate(const(char)[]) sink, string formatString) const @system
    {
        //sink(formatString);
        //formattedWrite(sink, formatString.replace("%r$", "%2$").replace("%c$", "%3$"), this.toRange!(Major.row), this.rlength, this.clength);
        foreach(i; 0 .. rlength)
            foreach(j; 0 .. clength)
                formattedWrite(sink, "%s, ", this[i, j]);
    }
    //*/
}


/**
要素がメモリ上に連続して存在するような行列を返します。

Example:
---
import extml.all;

void main(){
    auto mr = matrix!(2, 2, int, Major.row);    // 2x2, int型, 行優先
    assert(mr.array.equal([0, 0, 0, 0]));       // 初期化される

    mr[0, 1] = 1;
    mr[1, 0] = 2;
    mr[1, 1] = 3;
    assert(mr.array.equal([0, 1, 2, 3]));       // 行優先順


    auto mc = matrix!(2, 2, int, Major.column); // 2x2, int型, 列優先
    assert(mc.array.equal([0, 0, 0, 0]));       // 初期化される

    mc[0, 1] = 1;
    mc[1, 0] = 2;
    mc[1, 1] = 3;
    assert(mc.array.equal([0, 2, 1, 3]));       // 列優先順


    auto minit = matrix!(2, 2)(12);             // 2x2, int型, 全要素12
    assert(mc.array.equal([12, 12, 12, 12]));   // 全要素12で初期化されている

    minit[0, 0] = 0;
    minit[0, 1] = 1;
    minit[1, 0] = 2;
    minit[1, 1] = 3;

  version(DinearColumnMajor)
    assert(minit.array.equal([0, 2, 1, 3]));    // 列優先順
  else
    assert(minit.array.equal([0, 1, 2, 3]));    // 行優先順
}
---
*/
@property
auto matrix(size_t r, size_t c, T, Major mjr = defaultMajor)()
if(r != 0 && c != 0 && isScalar!T)
{
    static struct Matrix()
    {
        enum bool isColumnMajor = mjr == Major.column;
        enum bool isRowMajor    = mjr == Major.row;
        enum major = mjr;
        enum size_t rlength = r;
        enum size_t clength = c;


        ref inout(T) opIndex(size_t i, size_t j) inout
        in{
            assert(i < rlength);
            assert(j < clength);
        }
        body{
            static if(mjr == Major.row)
                return _array[i * clength + j];
            else
                return _array[j * rlength + i];
        }


        @property
        auto ref array() pure nothrow @safe inout
        {
            return _array[];
        }


        @property
        static typeof(this) init() pure nothrow @safe
        {
            return typeof(this)(new T[r * c]);
        }


        mixin ExpressionTemplateOperators!(ETOSpec.all);
    
      static if(is(typeof(_opIndex)))
        mixin(_opIndex);


      private:
        T[] _array;
    }

    return Matrix!().init;
}


///ditto
auto matrix(size_t r, size_t c, Major mjr = defaultMajor, T)(T init)
if(r != 0 && c != 0 && isScalar!T)
{
    auto m = matrix!(r, c, T, mjr)();
    m._array[] = init;
    return m;
}



auto matrix(size_t r, size_t c, Major mjr = defaultMajor, R)(R v)
if(r != 0 && c != 0 && isScalar!(Unqual!(std.range.ElementType!R)))
{
    auto m = matrix!(r, c, Unqual!(std.range.ElementType!R), mjr)();
    auto ar = m.array;
    ar.put(v);
    return m;
}


@property
auto rvector(size_t c, T, Major mjr = Major.row)()
if(c != 0 && isScalar!T)
{
    return matrix!(1, c, T, mjr);
}


auto rvector(size_t c, Major mjr = Major.row, T)(T init)
if(c != 0 && isScalar!T)
{
    return matrix!(1, c, mjr)(init);
}


auto rvector(size_t c, Major mjr = Major.row, R)(R v)
if(c != 0 && isScalar!(Unqual!(std.range.ElementType!R)))
{
    return matrix!(1, c, mjr)(v);
}


@property
auto cvector(size_t r, T, Major mjr = Major.column)()
if(r != 0 && isScalar!T)
{
    return matrix!(r, 1, T, mjr);
}


auto cvector(size_t r, Major mjr = Major.column, T)(T init)
if(r != 0 && isScalar!T)
{
    return matrix!(r, 1, mjr)(init);
}


auto cvector(size_t r, Major mjr = Major.column, R)(R v)
if(r != 0 && isScalar!(Unqual!(std.range.ElementType!R)))
{
    return matrix!(r, 1, mjr)(v);
}


unittest{
    scope(failure) writefln("Unittest failure : ", __FILE__, " : ", __LINE__);
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}


    auto mr = matrix!(2, 2, int, Major.row);    // 2x2, int型, 行優先
    assert(mr.array.equal([0, 0, 0, 0]));       // 初期化される

    mr[0, 1] = 1;
    mr[1, 0] = 2;
    mr[1, 1] = 3;
    assert(mr.array.equal([0, 1, 2, 3]));       // 行優先順


    auto mc = matrix!(2, 2, int, Major.column); // 2x2, int型, 列優先
    assert(mc.array.equal([0, 0, 0, 0]));       // 初期化される

    mc[0, 1] = 1;
    mc[1, 0] = 2;
    mc[1, 1] = 3;
    assert(mc.array.equal([0, 2, 1, 3]));       // 列優先順


    auto minit = matrix!(2, 2)(12);             // 2x2, int型, 全要素12
    assert(minit.array.equal([12, 12, 12, 12]));   // 全要素12で初期化されている

    minit[0, 0] = 0;
    minit[0, 1] = 1;
    minit[1, 0] = 2;
    minit[1, 1] = 3;

  version(DinearColumnMajor)
    assert(minit.array.equal([0, 2, 1, 3]));    // 列優先順
  else
    assert(minit.array.equal([0, 1, 2, 3]));    // 行優先順
}
unittest{
    auto m = matrix!(1, 1)(3);
    int a = m;
    assert(a == 3);
}
unittest{
    auto rv1 = rvector!(3)(3);
    static assert(rv1.rlength == 1);
    static assert(rv1.clength == 3);
    static assert(rv1.length  == 3);
    assert(equal(rv1.toFlatten, [3, 3, 3]));

    auto rv2 = rvector!(4)(rv1.toFlatten);
    static assert(rv2.rlength == 1);
    static assert(rv2.clength == 4);
    static assert(rv2.length  == 4);
    assert(equal(rv2.toFlatten, [3, 3, 3, int.init]));

    auto cv1 = cvector!(3)(3);
    static assert(cv1.rlength == 3);
    static assert(cv1.clength == 1);
    static assert(cv1.length  == 3);
    assert(equal(cv1.toFlatten, [3, 3, 3]));

    auto cv2 = cvector!(4)(cv1.toFlatten);
    static assert(cv2.rlength == 4);
    static assert(cv2.clength == 1);
    static assert(cv2.length  == 4);
    assert(equal(cv2.toFlatten, [3, 3, 3, int.init]));
}


auto toRange(Major mjr = defaultMajor, A)(A mat)
if(isMatrix!A)
{
  static if(mjr == Major.column)
    return mat.transpose.toRange!(Major.row);
  else
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
            size_t _cb = A.clength;
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
        size_t _rb = A.rlength;
    }


    return ToRange!()(mat);
  }
}

unittest{
    scope(failure) writefln("Unittest failure : ", __FILE__, " : ", __LINE__);
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}


    auto rm33 = matrix!(3, 3, int, Major.row);
    rm33[0, 0] = 1; rm33[0, 1] = 2; rm33[0, 2] = 3;
    assert(equal!"equal(a, b)"(rm33.toRange, [[1, 2, 3], [0, 0, 0], [0, 0, 0]]));
}



/**
行列をレンジにします
*/
@property
auto toFlatten(Major mjr = Major.row, A)(A mat)
if(isMatrix!A)
{
    alias ElementType!A E;

    static struct ToFlatten()
    {
        @property
        auto ref front()
        {
            static if(mjr == Major.row)
                return _mat[_f / A.clength, _f % A.clength];
            else
                return _mat[_f % A.rlength, _f / A.rlength];
        }


        @property
        auto ref back()
        {
            static if(mjr == Major.row)
                return _mat[_b / A.clength, _b % A.clength];
            else
                return _mat[_b % A.rlength, _b / A.rlength];
        }


        auto ref opIndex(size_t i)
        in{
            assert(_f + i < _b);
        }body{
            i += _f;
            static if(mjr == Major.row)
                return _mat[i / A.clength, i % A.clength];
            else
                return _mat[i % A.rlength, i / A.rlength];
        }


        static if(hasAssignableElements!A)
        {
            @property
            void front(E v)
            {
                static if(mjr == Major.row)
                    _mat[_f / A.clength, _f % A.clength] = v;
                else
                    _mat[_f % A.rlength, _f / A.rlength] = v;
            }


            @property
            void back(E v)
            {
                static if(mjr == Major.row)
                    _mat[_b / A.clength, _b % A.clength] = v;
                else
                    _mat[_b % A.rlength, _b / A.rlength] = v;
            }


            void opIndexAssign(E v, size_t i)
            in{
                assert(_f + i < _b);
            }
            body{
                i += _f;
                static if(mjr == Major.row)
                    _mat[i / A.clength, i % A.clength] = v;
                else
                    _mat[i % A.rlength, i / A.rlength] = v;
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
        size_t _f = 0, _b = A.rlength * A.clength;
    }

    return ToFlatten!()(mat);
}
unittest{
    scope(failure) writefln("Unittest failure : ", __FILE__, " : ", __LINE__);
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}


    auto rm33 = matrix!(3, 3, int, Major.row);
    rm33[0, 0] = 1; rm33[0, 1] = 2; rm33[0, 2] = 3;
    //writeln(rm33.array.length;)
    assert(rm33.array == [1, 2, 3, 0, 0, 0, 0, 0, 0]);

    alias Rt1 = typeof(rm33.toFlatten!(Major.row));
    static assert(isRandomAccessRange!(Rt1));
    static assert(std.range.hasLvalueElements!(Rt1));
    static assert(std.range.hasAssignableElements!(Rt1));
    static assert(hasLength!(Rt1));
    assert(equal(rm33.toFlatten!(Major.row), rm33.array));

    auto cm33 = matrix!(3, 3, int, Major.column);
    cm33[0, 0] = 1; cm33[0, 1] = 2; cm33[0, 2] = 3;
    assert(cm33.array == [1, 0, 0, 2, 0, 0, 3, 0, 0]);
    assert(equal(rm33.toFlatten!(Major.column), cm33.array));
}


/**
レンジから行列を作ります。
*/
auto toMatrix(size_t r, size_t c, Major mjr = defaultMajor, R)(R range)
if(isRandomAccessRange!R && isScalar!(Unqual!(std.range.ElementType!R)) && (mjr == Major.row ? (c != 0) : (r != 0)))
{
  //static if(mjr == Major.column)
    //return range.toMatrix!(c, r, Major.row).transpose;
  //else
  //{
    alias E = Unqual!(std.range.ElementType!R);

    static struct ToMatrix()
    {
        enum size_t rlength = r;
        enum size_t clength = c;


        auto ref opIndex(size_t i, size_t j) inout
        in{
            assert(i < rlength || rlength == 0);
            assert(j < clength || clength == 0);

          static if(hasLength!R && mjr == Major.row)
            assert(i * clength + j < this._range.length);

          static if(hasLength!R && mjr == Major.column)
            assert(j * rlength + i < this._range.length);
        }
        body{
          static if(mjr == Major.row)
            return this._range[i * clength + j];
          else
            return this._range[j * rlength + i];
        }


      static if(rlength == 0 || clength == 0)
      {
        static struct Instantiate(size_t r, size_t c)
        {
            enum bool isValid = isEqOrEitherEq0(ToMatrix.rlength, r)
                             && isEqOrEitherEq0(ToMatrix.clength, c);

          static if(isValid)
          {
            enum size_t rlength = r == 0 ? ToMatrix.rlength : r;
            enum size_t clength = c == 0 ? ToMatrix.clength : c;
          }
        }
      }


        mixin ExpressionTemplateOperators!(ETOSpec.all);

      static if(is(typeof(_opIndex)))
        mixin(_opIndex);

      private:
        R _range;
    }

    return ToMatrix!()(range);
  //}
}

unittest{
    scope(failure) writefln("Unittest failure : ", __FILE__, " : ", __LINE__);
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}


    auto r = iota(4);
    auto mr = toMatrix!(2, 2)(r);
    assert(mr.toFlatten.equal([0, 1, 2, 3]));

    auto mc = r.toMatrix!(2, 2, Major.column);
    assert(mc.toFlatten.equal([0, 2, 1, 3]));
}
unittest{
    scope(failure) writefln("Unittest failure : ", __FILE__, " : ", __LINE__);
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
    import std.stdio;
    assert(mc.toFlatten!(Major.column).equal([0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                                            0, 0, 0, 0, 0, 15]));
}
unittest{
    scope(failure) writefln("Unittest failure : ", __FILE__, " : ", __LINE__);
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}

    auto mem = [0, 1, 2, 3];
    auto r1 = mem.toMatrix!(0, 1, Major.row);
    assert(equal(r1.instantiate!(3, 1).toFlatten, mem[0 .. 3]));

    mem ~= 4;
    auto r2 = mem.toMatrix!(0, 2, Major.row);
    assert(r2[2, 0] == 4);

    auto c1 = mem.toMatrix!(1, 0, Major.column);
    assert(equal(instantiate!(1, 3)(c1).toFlatten, mem[0 .. 3]));
}



auto toMatrix(size_t r, size_t c, Major mjr = defaultMajor, R)(R range)
if(isRandomAccessRange!R && isRandomAccessRange!(Unqual!(std.range.ElementType!R)) && isScalar!(Unqual!(std.range.ElementType!(Unqual!(std.range.ElementType!R)))))
{
    static struct ToMatrix()
    {
        enum rlength = r;
        enum clength = c;


        auto ref opIndex(size_t i, size_t j) inout
        in{
            assert(i < rlength || rlength == 0);
            assert(j < clength || clength == 0);

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


      static if(rlength == 0 || clength == 0)
      {
        static struct Instantiate(size_t r, size_t c)
        {
            enum bool isValid = isEqOrEitherEq0(r, ToMatrix.rlength)
                             && isEqOrEitherEq0(c, ToMatrix.clength);

          static if(isValid)
          {
            enum size_t rlength = r == 0 ? ToMatrix.rlength : r;
            enum size_t clength = c == 0 ? ToMatrix.clength : c;
          }
        }
      }

        mixin ExpressionTemplateOperators!(ETOSpec.all);

      static if(is(typeof(_opIndex)))
        mixin(_opIndex);

      private:
        R _range;
    }
  
    return ToMatrix!()(range);
}
unittest{
    scope(failure) writefln("Unittest failure : ", __FILE__, " : ", __LINE__);
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}

    auto arr = [[0, 1], [2, 3], [4, 5]];
    auto r1 = arr.toMatrix!(3, 2, Major.row);
    assert(equal!"equal(a, b)"(r1.toRange!(Major.row), arr));

    auto r2 = arr.toMatrix!(1, 1, Major.row);
    assert(r2[0] == 0);

    auto r3 = arr.toMatrix!(0, 2, Major.row);
    assert(r3.instantiate!(3, 2) == r1);

    auto r4 = arr.toMatrix!(2, 0, Major.row);
    assert(equal(r4.instantiate!(2, 2).toFlatten!(Major.row), [0, 1, 2, 3]));

    auto r5 = arr.toMatrix!(0, 0, Major.row);
    assert(r5 == r1);
}
unittest{
    scope(failure) writefln("Unittest failure : ", __FILE__, " : ", __LINE__);
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}

    auto arr = [[0, 1], [2, 3], [4, 5]];
    auto r1 = arr.toMatrix!(2, 3, Major.column);
    assert(equal!"equal(a, b)"(r1.toRange!(Major.column), arr));
    assert(r1[0, 0] == 0); assert(r1[0, 1] == 2); assert(r1[0, 2] == 4);
    assert(r1[1, 0] == 1); assert(r1[1, 1] == 3); assert(r1[1, 2] == 5);

    auto r2 = arr.toMatrix!(1, 1, Major.column);
    assert(r2[0] == 0);

    auto r3 = arr.toMatrix!(0, 3, Major.column);
    assert(r3.instantiate!(2, 3) == r1);

    auto r4 = arr.toMatrix!(2, 0, Major.column);
    assert(equal(r4.instantiate!(2, 2).toFlatten!(Major.column), [0, 1, 2, 3]));

    auto r5 = arr.toMatrix!(0, 0, Major.column);
    assert(r5 == r1);
}


/**
インターフェース
*/
interface IMatrix(T, size_t r, size_t c)
{
    enum rlength = r;
    enum clength = c;

  static if(rlength == 0 || clength == 0)
  {
    static struct Instantiate(size_t i, size_t j)
    {
        enum isValid = (IMatrix.rlength == 0 ? true : IMatrix.rlength == i)
                    && (IMatrix.clength == 0 ? true : IMatrix.clength == i);

        static if(isValid)
        {
            enum size_t rlength = IMatrix.rlength;
            enum size_t clength = IMatrix.clength;
        }
    }
  }

    T opIndex(size_t i, size_t j)
    in{
        assert(i < rlength || rlength == 0);
        assert(j < clength || clength == 0);
    }


    T opIndex(size_t i, size_t j) const 
    in{
        assert(i < rlength || rlength == 0);
        assert(j < clength || clength == 0);
    }


    T opIndex(size_t i, size_t j) immutable
    in{
        assert(i < rlength || rlength == 0);
        assert(j < clength || clength == 0);
    }
}



/**
参照型オブジェクトにして返します。
*/
@property
IMatrix!(ElementType!A, A.rlength, A.clength)
        toInterface(A)(A mat)
if(isMatrix!A)
{
    alias ElementType!A E;

    static class MatrixObj() : IMatrix!(ElementType!A, A.rlength, A.clength)
    {
        this(A mat)
        {
            _mat = mat;
        }


        T opIndex(size_t i, size_t j)
        in{
            assert(i < rlength);
            assert(j < clength);
        }
        body{
            return _mat[i, j];
        }


        T opIndex(size_t i, size_t j) const
        in{
            assert(i < rlength);
            assert(j < clength);
        }
        body{
            return _mat[i, j];
        }


        T opIndex(size_t i, size_t j) immutable
        in{
            assert(i < rlength);
            assert(j < clength);
        }
        body{
            return _mat[i, j];
        }


      private:
        A _mat;
    }


    return new MatrixObj!()(mat);
}




/**
単位行列を作ります。

Example:
---
void main(){
    auto id = identity!int;

    static assert(typeof(id).Instantiate!(4, 4).isValid);
    static assert(!typeof(id).Instantiate!(1, 3).isValid);

    auto m1 = Matrix2i.init;
    auto r1 = m1.toRange;
    r1.put([0, 1, 2, 3]);
    static assert(isRandomAccessRange!(typeof(r1)));
    static assert(std.range.hasLvalueElements!(typeof(r1)));
    assert(equal((m1 * id).range, [0, 1, 2, 3]));

    auto id2 = id + id;
    static assert(isMatrix!(typeof(id2)));
    static assert(typeof(id2).Instantiate!(4, 4).isValid);
    static assert(isFlexibleMatrix!(typeof(id2)));


    static assert(isFlexibleMatrix!(typeof(id.instantiate!(0, 2))));
    auto id3 = id.instantiate!(0, 2) * id;
    static assert(id3.rlength == 2);
    static assert(id3.clength == 2);
    static assert(!isFlexibleMatrix!(typeof(id3)));
    assert(equal(id3.toRange, [1, 0, 0, 1]));

    auto ins = id2.instantiate!(2, 2);
    static assert(isMatrix!(typeof(ins)));
    static assert(!isFlexibleMatrix!(typeof(ins)));
    assert(equal(ins.toRange, [2, 0, 0, 2]));
}
---
*/
@property
auto identity(E)()if(isScalar!E)
{
    static struct Identity()
    {
        enum rlength = 0;
        enum clength = 0;


        E opIndex(size_t i, size_t j) inout
        {
            return (i == j) ? (cast(E)1) : (cast(E)0);
        }


      static struct Instantiate(size_t i, size_t j)
      {
        static if(i == 0 && j != 0)
        {
            enum isValid = true;
            enum rlength = j;
            enum clength = j;
        }
        else static if(j == 0 && i != 0)
        {
            enum isValid = true;
            enum rlength = i;
            enum clength = i;
        }
        else
        {
            enum isValid = i == j;
            static if(isValid)
            {
                enum rlength = i;
                enum clength = j;
            }
        }
      }

        mixin ExpressionTemplateOperators!(ETOSpec.all);

      static if(is(typeof(_opIndex)))
        mixin(_opIndex);
    }

    static assert(isMatrix!(Identity!()));

    return Identity!()();
}


unittest{
    scope(failure) writefln("Unittest failure : ", __FILE__, " : ", __LINE__);
    scope(success) {writefln("Unittest success :%s(%s)", __FILE__, __LINE__); stdout.flush();}

    auto id = identity!int;

    static assert(typeof(id).Instantiate!(4, 4).isValid);
    static assert(!typeof(id).Instantiate!(1, 3).isValid);

    auto m1 = Matrix2i.init;
    auto r1 = m1.toFlatten!(Major.row);
    r1.put([0, 1, 2, 3]);
    static assert(std.range.hasLvalueElements!(typeof(r1)));
    assert(equal((m1 * id).toFlatten!(Major.row), [0, 1, 2, 3]));

    auto id2 = id + id;
    static assert(isMatrix!(typeof(id2)));
    static assert(typeof(id2).Instantiate!(4, 4).isValid);
    static assert(isFlexibleMatrix!(typeof(id2)));


    static assert(isFlexibleMatrix!(typeof(id.instantiate!(0, 2))));
    auto id3 = id.instantiate!(0, 2) * id;
    static assert(id3.rlength == 2);
    static assert(id3.clength == 2);
    static assert(!isFlexibleMatrix!(typeof(id3)));
    assert(equal(id3.toFlatten, [1, 0, 0, 1]));

    auto ins = id2.instantiate!(2, 2);
    static assert(isMatrix!(typeof(ins)));
    static assert(!isFlexibleMatrix!(typeof(ins)));
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
        enum rlength = 0;
        enum clength = 0;


        E opIndex(size_t i, size_t j) inout
        {
            return cast(E)1;
        }


        static struct Instantiate(size_t i, size_t j)
        {
            enum isValid = true;
            enum rlength = i;
            enum clength = j;
        }


        mixin ExpressionTemplateOperators!(ETOSpec.all & ~ETOSpec.matrixMulScalar & ~ETOSpec.scalarMulMatrix);

        static if(is(typeof(_opIndex)))
        mixin(_opIndex);


        auto opBinary(string op : "*", S)(S s)
        if(isScalar!S)
        {
            static assert(isValidOperator!(typeof(this), op, S));
            return ns(s);
        }


        auto opBinary(string op : "*", S)(const S s)
        if(isScalar!S)
        {
            static assert(isValidOperator!(typeof(this), op, const(S)));
            return ns(s);
        }


        auto opBinary(string op : "*", S)(immutable S s)
        if(isScalar!S)
        {
            static assert(isValidOperator!(typeof(this), op, immutable(S)));
            return ns(s);
        }


        auto opBinaryRight(string op : "*", S)(S s)
        if(isScalar!S)
        {
            static assert(isValidOperator!(S, op, typeof(this)));
            return ns(s);
        }


        auto opBinaryRight(string op : "*", S)(const S s)
        if(isScalar!S)
        {
            static assert(isValidOperator!(typeof(this), op, const(S)));
            return ns(s);
        }


        auto opBinaryRight(string op : "*", S)(immutable S s)
        if(isScalar!S)
        {
            static assert(isValidOperator!(typeof(this), op, immutable(S)));
            return ns(s);
        }


        static struct Ns(E)
        {
            enum rlength = 0;
            enum clength = 0;

            E opIndex(size_t i, size_t j) inout
            {
                return e;
            }


            static struct Instantiate(size_t i, size_t j)
            {
                enum isValid = true;
                enum rlength = i;
                enum clength = j;
            }


            mixin ExpressionTemplateOperators!(ETOSpec.all & ~ETOSpec.matrixMulScalar & ~ETOSpec.scalarMulMatrix);

            static if(is(typeof(_opIndex)))
            mixin(_opIndex);


            auto opBinary(string op : "*", S)(S s)
            if(isScalar!S)
            {
                static assert(isValidOperator!(typeof(this), op, const(S)));
                return ns(e * s);
            }


            auto opBinary(string op : "*", S)(const S s)
            if(isScalar!S)
            {
                static assert(isValidOperator!(typeof(this), op, const(S)));
                return ns(e * s);
            }


            auto opBinary(string op : "*", S)(immutable S s)
            if(isScalar!S)
            {
                static assert(isValidOperator!(typeof(this), op, immutable(S)));
                return ns(e * s);
            }


            auto opBinaryRight(string op : "*", S)(S s)
            if(isScalar!S)
            {
                static assert(isValidOperator!(S, op, typeof(this)));
                return ns(s * e);
            }


            auto opBinaryRight(string op : "*", S)(const S s)
            if(isScalar!S)
            {
                static assert(isValidOperator!(typeof(this), op, const(S)));
                return ns(s * e);
            }


            auto opBinaryRight(string op : "*", S)(immutable S s)
            if(isScalar!S)
            {
                static assert(isValidOperator!(typeof(this), op, immutable(S)));
                return ns(s * e);
            }

          private:
            E e;
        }

      private:
        static auto ns(S)(S s)
        {
            return Ns!S(s);
        }
    }

    static assert(isMatrix!(Ones!()));

    return Ones!().init;
}
unittest{
    auto m1 = ones!float;
    assert(m1[0, 1] == 1);

    auto m3 = m1 * 3;
    assert(m3[0, 1] == 3);

    auto m9 = m3 * 3;
    assert(m9[0, 1] == 9);
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
        enum size_t rlength = A.clength;
        enum size_t clength = A.rlength;


        auto ref opIndex(size_t i, size_t j)
        in{
            assert(i < rlength);
            assert(j < clength);
        }
        body{
            return mat[j, i];
        }


        auto ref opIndex(size_t i, size_t j) const
        in{
            assert(i < rlength);
            assert(j < clength);
        }
        body{
            return mat[j, i];
        }


        auto ref opIndex(size_t i, size_t j) immutable
        in{
            assert(i < rlength);
            assert(j < clength);
        }
        body{
            return mat[j, i];
        }


      static if(extml.core.hasAssignableElements!A)
      {
        void opIndexAssign(E)(E e, size_t i, size_t j) @system
        in{
            assert(i < rlength);
            assert(j < clength);
        }
        body{
            pragma(msg, __FILE__, "(", __LINE__, "): @system...");
            mat[j, i] = e;
        }
      }


      static if(rlength == 0 || clength == 0)
      {
        template Instantiate(size_t r, size_t c)
        {
            alias Instantiate = A.Instantiate!(c, r);
        }
      }


        mixin ExpressionTemplateOperators!(ETOSpec.all);

      static if(is(typeof(_opIndex)))
        mixin(_opIndex);

      private:
        A mat;
    }

    static assert(isMatrix!(Transposed!()));

    return Transposed!()(mat);
}
unittest{
    auto m = matrix!(2, 2, int)();
    m[0, 0] = 0; m[0, 1] = 1;
    m[1, 0] = 2; m[1, 1] = 3;

    auto t = transpose(m);
    static assert(isMatrix!(typeof(t)));
    assert(equal(t.toFlatten!(Major.row), [0, 2, 1, 3]));

    t[0, 1] = 100;
    assert(m[1, 0] == 100);
}


/**
部分行列を返します
*/
auto sub(alias rArray, alias cArray, A)(A mat)
if(isArray!(typeof(rArray)) && isArray!(typeof(cArray)) && isMatrix!A)
{
    static assert(rArray.length != 0 || (A.rlength == 0 && isFlexibleMatrix!A));
    static assert(cArray.length != 0 || (A.clength == 0 && isFlexibleMatrix!A));

    static assert(A.rlength == 0 || rArray.find!"a>=b"(A.rlength).empty);
    static assert(A.clength == 0 || cArray.find!"a>=b"(A.clength).empty);

  static if(rArray.length == 0 && cArray.length == 0)
    return mat;
  else
  {
    static struct Sub()
    {
        enum rlength = rArray.length;
        enum clength = cArray.length;


        auto ref opIndex(size_t i, size_t j)
        in{
            assert(i < rlength);
            assert(j < clength);
        }
        body{
            static if(rlength != 0 && clength != 0)
            return _mat[rArray[i], cArray[j]];
          else static if(rlength != 0)
            return _mat[rArray[i], j];
          else
            return _mat[i, cArray[j]];
        }


        auto ref opIndex(size_t i, size_t j) const
        in{
            assert(i < rlength || rlength == 0);
            assert(j < clength || clength == 0);
        }
        body{
          static if(rlength != 0 && clength != 0)
            return _mat[rArray[i], cArray[j]];
          else static if(rlength != 0)
            return _mat[rArray[i], j];
          else
            return _mat[i, cArray[j]];
        }


        auto ref opIndex(size_t i, size_t j) immutable
        in{
            assert(i < rlength);
            assert(j < clength);
        }
        body{
            static if(rlength != 0 && clength != 0)
            return _mat[rArray[i], cArray[j]];
          else static if(rlength != 0)
            return _mat[rArray[i], j];
          else
            return _mat[i, cArray[j]];
        }


        mixin ExpressionTemplateOperators!(ETOSpec.all);
        static if(is(typeof(_opIndex)))
            mixin(_opIndex);


        static if(rlength == 0 || clength == 0)
        static struct Instantiate(size_t i, size_t j)
        {
            enum isValid = isEqOrEitherEq0(i, Sub.rlength)
                        && isEqOrEitherEq0(j, Sub.clength);

            static if(isValid)
            {
                enum rlength = i == 0 ? Sub.rlength : i;
                enum clength = j == 0 ? Sub.clength : j;
            }
        }


      private:
        A _mat;
    }

    return Sub!()(mat);
  }
}
unittest{
    auto m1 = [[0, 1], [2, 3], [4, 5]].toMatrix!(3, 2, Major.row);
    auto s1 = m1.sub!([1, 2], [0]);
    static assert(s1.rlength == 2);
    static assert(s1.clength == 1);

    assert(s1[0, 0] == 2);
    assert(s1[1, 0] == 4);


    auto m2 = [[0, 1], [2, 3], [4, 5]].toMatrix!(0, 0, Major.row);
    auto s2 = sub!((size_t[]).init, (size_t[]).init)(m2);
    assert(m1 == s2.instantiate!(3, 2));


    auto m3 = sub!([0, 0, 0], [0, 0])(identity!float);
    assert(m3 == ones!float);


    auto m4 = sub!([0], (size_t[]).init)(identity!float);
    static assert(m4.Instantiate!(1, 3).isValid);
    static assert(!m4.Instantiate!(3, 3).isValid);
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

    return sub!(arrGen(r, A.rlength), arrGen(c, A.clength))(mat);
}
unittest{
    auto a = swizzle!("abab", "baba")(identity!int);
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
@property
ElementType!A trace(A)(A mat)
if(isMatrix!A && (!isFlexibleMatrix!A && A.rlength == A.clength))
{
    alias ElementType!A T;
    T sum = cast(T)0;
    foreach(i; 0 .. A.rlength)
        sum += mat[i, i];
    return sum;
}
unittest{
    auto m = matrix!(2, 2, int)();
    m[0, 0] = 0; m[0, 1] = 1;
    m[1, 0] = 2; m[1, 1] = 3;

    auto tr = m.trace;
    assert(tr == 3);
}



///ditto
@property
ElementType!A trace(size_t rc, A)(A mat)
if(isMatrix!A && isFlexibleMatrix!A)
{
    static assert(A.Instantiate!(rc, rc).isValid);
    return mat.instantiate!(rc, rc).trace;
}


auto dot(V1, V2)(V1 vec1, V2 vec2)
if(isVector!V1 && isVector!V2)
{
    static if(V1.rlength == 1)
    {
        static if(V2.clength == 1)
            return (vec1 * vec2)[0];
        else
            return (vec1 * vec2.transpose)[0];
    }
    else
    {
        static if(V2.rlength == 1)
            return (vec2 * vec1)[0];    //vec1.transpose * vec2.tranpose == vec2 * vec1
        else
            return (vec1.transpose * vec2)[0];
    }
}
unittest{
    auto rv = rvector!3([0, 1, 2][]),
         cv = cvector!3([1, 2, 3][]);

    assert((rv * cv)[0] == 8);  //8
    assert(rv.dot(cv) == 8);
    assert(cv.dot(rv) == 8);

    assert(rv.dot(rv) == 5);
    assert(cv.dot(cv) == 14);
}


/*
auto cross(Major mjr = defaultMajor, V1, V2)(V1 vec1, V2 vec2)
*/

/**
直積
*/
auto cartesian(V1, V2)(V1 vec1, V2 vec2)
if(isVector!V1 && isVector!V2)
{
    static struct Cartesian()
    {
        enum rlength = V1.length;
        enum clength = V2.length;


        auto opIndex(size_t i, size_t j){ return _vec1[i] * _vec2[j]; }
        auto opIndex(size_t i, size_t j) const { return _vec1[i] * _vec2[j]; }
        auto opIndex(size_t i, size_t j) immutable { return _vec1[i] * _vec2[j]; }

      static if(rlength == 0 || clength == 0)
        static struct Instantiate(size_t i, size_t j)
        {
            enum isValid = isEqOrEitherEq0(i, rlength)
                        && isEqOrEitherEq0(j, clength);

            static if(isValid)
            {
                enum rlength = rlength == 0 ? r : rlength;
                enum clength = clength == 0 ? c : clength;
            }
        }


        static if(is(typeof(_opIndex)))
        mixin(_opIndex);

      private:
        V1 _vec1;
        V2 _vec2;
    }

    return Cartesian!()(vec1, vec2);
}
unittest{
    auto v1 = [0, 1, 2, 3].toMatrix!(3, 1);
    auto v2 = [2, 3, 4, 5].toMatrix!(1, 2);

    assert(v1.cartesian(v2) == v1 * v2);
}


/**
matrix!(...)で生成される型です。
フォーマットは、
Matrix(R)x(C)(Type)(Major)
となっており、
1 <= R <= 4
1 <= C <= 4

Type : int -> "i"
       uint -> "ui"
       ...
       float -> "f"
       cfloat -> "cf"
       ...
Major : "" <- デフォルト
        "r"　<- 行優先
        "c" <- 列優先

さらにR == Cのときは
Matrix(R)(Type)(Major)

R == 1のときは行ベクトルなので
Rvector(C)(Type)(Major)

C == 1のときは列ベクトルなので
Cvector(R)(Type)(Major)

とも定義されます。
*/
/+
import std.stdio;

void main(){
  writeln(genTypes);
}

auto genTypes(){
    import std.string;

    string[][] range(string[string] hstr)
    {
        typeof(return) dst;
        foreach(k, v; hstr)
            dst ~= [k, v];
        return dst;
    }


    string str;

    auto sf = ["byte" : "b",
               "ubyte" : "ub",
               "short" : "s",
               "ushort" : "us",
               "int" : "i",
               "uint" : "ui",
               "long" : "l",
               "ulong" : "ul",
               "float" : "f",
               "double" : "d",
               "real" : "r",
               "cfloat" : "cf",
               "cdouble" : "cd",
               "creal" : "cr"];

    
    auto dim = [[1, 1], [1, 2], [1, 3], [1, 4],
                [2, 1], [2, 2], [2, 3], [2, 4],
                [3, 1], [3, 2], [3, 3], [3, 4],
                [4, 1], [4, 2], [4, 3], [4, 4]];

    auto kv = range(sf);
    auto rc = [ ["defaultMajor", ""],
                ["Major.row", "r"],
                ["Major.column", "c"]
                ];

    foreach(d; dim)
        foreach(k; kv)
            foreach(r; rc)
            {
                if(d[0] == 1 && d[1] == 1)
                    continue;

                str ~= xformat(q{alias Matrix%1$sx%2$s%3$s%4$s = typeof(matrix!(%1$s, %2$s, %5$s, %6$s));} ~ "\n",
                                      d[0], d[1], k[1], r[1], k[0], r[0]);

                if(d[0] == d[1]){
                    str ~= xformat(q{alias Matrix%1$s%3$s%4$s = typeof(matrix!(%1$s, %2$s, %5$s, %6$s));} ~ "\n",
                                      d[0], d[1], k[1], r[1], k[0], r[0]);
                }

                if(r[1] == "c" && d[1] == 1){
                    str ~= xformat(q{alias Cvector%1$s%2$s = typeof(matrix!(%1$s, 1, %3$s, %4$s));} ~ "\n",
                                           d[0], k[1], k[0], r[0]);
                }

                if(r[1] == "r" && d[0] == 1){
                    str ~= xformat(q{alias Rvector%1$s%2$s = typeof(matrix!(1, %1$s, %3$s, %4$s));} ~ "\n",
                                           d[1], k[1], k[0], r[0]);
                }
            }
    return str;
}

mixin(genTypes());
+/
template Matrix(T, size_t r, size_t c, Major mjr = defaultMajor)
{
    alias Matrix = typeof(matrix!(r, c, T, mjr));
}


template Rvector(T, size_t c)
{
    alias Rvector = typeof(matrix!(T, 1, c, Major.row));
}


template Cvector(T, size_t r)
{
    alias Cvector = typeof(matrix!(T, r, 1, Major.column));
}


alias Matrix1x1ui = typeof(matrix!(1, 1, uint, defaultMajor));
alias Matrix1ui = typeof(matrix!(1, 1, uint, defaultMajor));
alias Matrix1x1uir = typeof(matrix!(1, 1, uint, Major.row));
alias Matrix1uir = typeof(matrix!(1, 1, uint, Major.row));
alias Rvector1ui = typeof(matrix!(1, 1, uint, Major.row));
alias Matrix1x1uic = typeof(matrix!(1, 1, uint, Major.column));
alias Matrix1uic = typeof(matrix!(1, 1, uint, Major.column));
alias Cvector1ui = typeof(matrix!(1, 1, uint, Major.column));
alias Matrix1x1f = typeof(matrix!(1, 1, float, defaultMajor));
alias Matrix1f = typeof(matrix!(1, 1, float, defaultMajor));
alias Matrix1x1fr = typeof(matrix!(1, 1, float, Major.row));
alias Matrix1fr = typeof(matrix!(1, 1, float, Major.row));
alias Rvector1f = typeof(matrix!(1, 1, float, Major.row));
alias Matrix1x1fc = typeof(matrix!(1, 1, float, Major.column));
alias Matrix1fc = typeof(matrix!(1, 1, float, Major.column));
alias Cvector1f = typeof(matrix!(1, 1, float, Major.column));
alias Matrix1x1s = typeof(matrix!(1, 1, short, defaultMajor));
alias Matrix1s = typeof(matrix!(1, 1, short, defaultMajor));
alias Matrix1x1sr = typeof(matrix!(1, 1, short, Major.row));
alias Matrix1sr = typeof(matrix!(1, 1, short, Major.row));
alias Rvector1s = typeof(matrix!(1, 1, short, Major.row));
alias Matrix1x1sc = typeof(matrix!(1, 1, short, Major.column));
alias Matrix1sc = typeof(matrix!(1, 1, short, Major.column));
alias Cvector1s = typeof(matrix!(1, 1, short, Major.column));
alias Matrix1x1d = typeof(matrix!(1, 1, double, defaultMajor));
alias Matrix1d = typeof(matrix!(1, 1, double, defaultMajor));
alias Matrix1x1dr = typeof(matrix!(1, 1, double, Major.row));
alias Matrix1dr = typeof(matrix!(1, 1, double, Major.row));
alias Rvector1d = typeof(matrix!(1, 1, double, Major.row));
alias Matrix1x1dc = typeof(matrix!(1, 1, double, Major.column));
alias Matrix1dc = typeof(matrix!(1, 1, double, Major.column));
alias Cvector1d = typeof(matrix!(1, 1, double, Major.column));
alias Matrix1x1cf = typeof(matrix!(1, 1, cfloat, defaultMajor));
alias Matrix1cf = typeof(matrix!(1, 1, cfloat, defaultMajor));
alias Matrix1x1cfr = typeof(matrix!(1, 1, cfloat, Major.row));
alias Matrix1cfr = typeof(matrix!(1, 1, cfloat, Major.row));
alias Rvector1cf = typeof(matrix!(1, 1, cfloat, Major.row));
alias Matrix1x1cfc = typeof(matrix!(1, 1, cfloat, Major.column));
alias Matrix1cfc = typeof(matrix!(1, 1, cfloat, Major.column));
alias Cvector1cf = typeof(matrix!(1, 1, cfloat, Major.column));
alias Matrix1x1ub = typeof(matrix!(1, 1, ubyte, defaultMajor));
alias Matrix1ub = typeof(matrix!(1, 1, ubyte, defaultMajor));
alias Matrix1x1ubr = typeof(matrix!(1, 1, ubyte, Major.row));
alias Matrix1ubr = typeof(matrix!(1, 1, ubyte, Major.row));
alias Rvector1ub = typeof(matrix!(1, 1, ubyte, Major.row));
alias Matrix1x1ubc = typeof(matrix!(1, 1, ubyte, Major.column));
alias Matrix1ubc = typeof(matrix!(1, 1, ubyte, Major.column));
alias Cvector1ub = typeof(matrix!(1, 1, ubyte, Major.column));
alias Matrix1x1b = typeof(matrix!(1, 1, byte, defaultMajor));
alias Matrix1b = typeof(matrix!(1, 1, byte, defaultMajor));
alias Matrix1x1br = typeof(matrix!(1, 1, byte, Major.row));
alias Matrix1br = typeof(matrix!(1, 1, byte, Major.row));
alias Rvector1b = typeof(matrix!(1, 1, byte, Major.row));
alias Matrix1x1bc = typeof(matrix!(1, 1, byte, Major.column));
alias Matrix1bc = typeof(matrix!(1, 1, byte, Major.column));
alias Cvector1b = typeof(matrix!(1, 1, byte, Major.column));
alias Matrix1x1cr = typeof(matrix!(1, 1, creal, defaultMajor));
alias Matrix1cr = typeof(matrix!(1, 1, creal, defaultMajor));
alias Matrix1x1crr = typeof(matrix!(1, 1, creal, Major.row));
alias Matrix1crr = typeof(matrix!(1, 1, creal, Major.row));
alias Rvector1cr = typeof(matrix!(1, 1, creal, Major.row));
alias Matrix1x1crc = typeof(matrix!(1, 1, creal, Major.column));
alias Matrix1crc = typeof(matrix!(1, 1, creal, Major.column));
alias Cvector1cr = typeof(matrix!(1, 1, creal, Major.column));
alias Matrix1x1i = typeof(matrix!(1, 1, int, defaultMajor));
alias Matrix1i = typeof(matrix!(1, 1, int, defaultMajor));
alias Matrix1x1ir = typeof(matrix!(1, 1, int, Major.row));
alias Matrix1ir = typeof(matrix!(1, 1, int, Major.row));
alias Rvector1i = typeof(matrix!(1, 1, int, Major.row));
alias Matrix1x1ic = typeof(matrix!(1, 1, int, Major.column));
alias Matrix1ic = typeof(matrix!(1, 1, int, Major.column));
alias Cvector1i = typeof(matrix!(1, 1, int, Major.column));
alias Matrix1x1ul = typeof(matrix!(1, 1, ulong, defaultMajor));
alias Matrix1ul = typeof(matrix!(1, 1, ulong, defaultMajor));
alias Matrix1x1ulr = typeof(matrix!(1, 1, ulong, Major.row));
alias Matrix1ulr = typeof(matrix!(1, 1, ulong, Major.row));
alias Rvector1ul = typeof(matrix!(1, 1, ulong, Major.row));
alias Matrix1x1ulc = typeof(matrix!(1, 1, ulong, Major.column));
alias Matrix1ulc = typeof(matrix!(1, 1, ulong, Major.column));
alias Cvector1ul = typeof(matrix!(1, 1, ulong, Major.column));
alias Matrix1x1l = typeof(matrix!(1, 1, long, defaultMajor));
alias Matrix1l = typeof(matrix!(1, 1, long, defaultMajor));
alias Matrix1x1lr = typeof(matrix!(1, 1, long, Major.row));
alias Matrix1lr = typeof(matrix!(1, 1, long, Major.row));
alias Rvector1l = typeof(matrix!(1, 1, long, Major.row));
alias Matrix1x1lc = typeof(matrix!(1, 1, long, Major.column));
alias Matrix1lc = typeof(matrix!(1, 1, long, Major.column));
alias Cvector1l = typeof(matrix!(1, 1, long, Major.column));
alias Matrix1x1us = typeof(matrix!(1, 1, ushort, defaultMajor));
alias Matrix1us = typeof(matrix!(1, 1, ushort, defaultMajor));
alias Matrix1x1usr = typeof(matrix!(1, 1, ushort, Major.row));
alias Matrix1usr = typeof(matrix!(1, 1, ushort, Major.row));
alias Rvector1us = typeof(matrix!(1, 1, ushort, Major.row));
alias Matrix1x1usc = typeof(matrix!(1, 1, ushort, Major.column));
alias Matrix1usc = typeof(matrix!(1, 1, ushort, Major.column));
alias Cvector1us = typeof(matrix!(1, 1, ushort, Major.column));
alias Matrix1x1r = typeof(matrix!(1, 1, real, defaultMajor));
alias Matrix1r = typeof(matrix!(1, 1, real, defaultMajor));
alias Matrix1x1rr = typeof(matrix!(1, 1, real, Major.row));
alias Matrix1rr = typeof(matrix!(1, 1, real, Major.row));
alias Rvector1r = typeof(matrix!(1, 1, real, Major.row));
alias Matrix1x1rc = typeof(matrix!(1, 1, real, Major.column));
alias Matrix1rc = typeof(matrix!(1, 1, real, Major.column));
alias Cvector1r = typeof(matrix!(1, 1, real, Major.column));
alias Matrix1x1cd = typeof(matrix!(1, 1, cdouble, defaultMajor));
alias Matrix1cd = typeof(matrix!(1, 1, cdouble, defaultMajor));
alias Matrix1x1cdr = typeof(matrix!(1, 1, cdouble, Major.row));
alias Matrix1cdr = typeof(matrix!(1, 1, cdouble, Major.row));
alias Rvector1cd = typeof(matrix!(1, 1, cdouble, Major.row));
alias Matrix1x1cdc = typeof(matrix!(1, 1, cdouble, Major.column));
alias Matrix1cdc = typeof(matrix!(1, 1, cdouble, Major.column));
alias Cvector1cd = typeof(matrix!(1, 1, cdouble, Major.column));
alias Matrix1x2ui = typeof(matrix!(1, 2, uint, defaultMajor));
alias Matrix1x2uir = typeof(matrix!(1, 2, uint, Major.row));
alias Rvector2ui = typeof(matrix!(1, 2, uint, Major.row));
alias Matrix1x2uic = typeof(matrix!(1, 2, uint, Major.column));
alias Matrix1x2f = typeof(matrix!(1, 2, float, defaultMajor));
alias Matrix1x2fr = typeof(matrix!(1, 2, float, Major.row));
alias Rvector2f = typeof(matrix!(1, 2, float, Major.row));
alias Matrix1x2fc = typeof(matrix!(1, 2, float, Major.column));
alias Matrix1x2s = typeof(matrix!(1, 2, short, defaultMajor));
alias Matrix1x2sr = typeof(matrix!(1, 2, short, Major.row));
alias Rvector2s = typeof(matrix!(1, 2, short, Major.row));
alias Matrix1x2sc = typeof(matrix!(1, 2, short, Major.column));
alias Matrix1x2d = typeof(matrix!(1, 2, double, defaultMajor));
alias Matrix1x2dr = typeof(matrix!(1, 2, double, Major.row));
alias Rvector2d = typeof(matrix!(1, 2, double, Major.row));
alias Matrix1x2dc = typeof(matrix!(1, 2, double, Major.column));
alias Matrix1x2cf = typeof(matrix!(1, 2, cfloat, defaultMajor));
alias Matrix1x2cfr = typeof(matrix!(1, 2, cfloat, Major.row));
alias Rvector2cf = typeof(matrix!(1, 2, cfloat, Major.row));
alias Matrix1x2cfc = typeof(matrix!(1, 2, cfloat, Major.column));
alias Matrix1x2ub = typeof(matrix!(1, 2, ubyte, defaultMajor));
alias Matrix1x2ubr = typeof(matrix!(1, 2, ubyte, Major.row));
alias Rvector2ub = typeof(matrix!(1, 2, ubyte, Major.row));
alias Matrix1x2ubc = typeof(matrix!(1, 2, ubyte, Major.column));
alias Matrix1x2b = typeof(matrix!(1, 2, byte, defaultMajor));
alias Matrix1x2br = typeof(matrix!(1, 2, byte, Major.row));
alias Rvector2b = typeof(matrix!(1, 2, byte, Major.row));
alias Matrix1x2bc = typeof(matrix!(1, 2, byte, Major.column));
alias Matrix1x2cr = typeof(matrix!(1, 2, creal, defaultMajor));
alias Matrix1x2crr = typeof(matrix!(1, 2, creal, Major.row));
alias Rvector2cr = typeof(matrix!(1, 2, creal, Major.row));
alias Matrix1x2crc = typeof(matrix!(1, 2, creal, Major.column));
alias Matrix1x2i = typeof(matrix!(1, 2, int, defaultMajor));
alias Matrix1x2ir = typeof(matrix!(1, 2, int, Major.row));
alias Rvector2i = typeof(matrix!(1, 2, int, Major.row));
alias Matrix1x2ic = typeof(matrix!(1, 2, int, Major.column));
alias Matrix1x2ul = typeof(matrix!(1, 2, ulong, defaultMajor));
alias Matrix1x2ulr = typeof(matrix!(1, 2, ulong, Major.row));
alias Rvector2ul = typeof(matrix!(1, 2, ulong, Major.row));
alias Matrix1x2ulc = typeof(matrix!(1, 2, ulong, Major.column));
alias Matrix1x2l = typeof(matrix!(1, 2, long, defaultMajor));
alias Matrix1x2lr = typeof(matrix!(1, 2, long, Major.row));
alias Rvector2l = typeof(matrix!(1, 2, long, Major.row));
alias Matrix1x2lc = typeof(matrix!(1, 2, long, Major.column));
alias Matrix1x2us = typeof(matrix!(1, 2, ushort, defaultMajor));
alias Matrix1x2usr = typeof(matrix!(1, 2, ushort, Major.row));
alias Rvector2us = typeof(matrix!(1, 2, ushort, Major.row));
alias Matrix1x2usc = typeof(matrix!(1, 2, ushort, Major.column));
alias Matrix1x2r = typeof(matrix!(1, 2, real, defaultMajor));
alias Matrix1x2rr = typeof(matrix!(1, 2, real, Major.row));
alias Rvector2r = typeof(matrix!(1, 2, real, Major.row));
alias Matrix1x2rc = typeof(matrix!(1, 2, real, Major.column));
alias Matrix1x2cd = typeof(matrix!(1, 2, cdouble, defaultMajor));
alias Matrix1x2cdr = typeof(matrix!(1, 2, cdouble, Major.row));
alias Rvector2cd = typeof(matrix!(1, 2, cdouble, Major.row));
alias Matrix1x2cdc = typeof(matrix!(1, 2, cdouble, Major.column));
alias Matrix1x3ui = typeof(matrix!(1, 3, uint, defaultMajor));
alias Matrix1x3uir = typeof(matrix!(1, 3, uint, Major.row));
alias Rvector3ui = typeof(matrix!(1, 3, uint, Major.row));
alias Matrix1x3uic = typeof(matrix!(1, 3, uint, Major.column));
alias Matrix1x3f = typeof(matrix!(1, 3, float, defaultMajor));
alias Matrix1x3fr = typeof(matrix!(1, 3, float, Major.row));
alias Rvector3f = typeof(matrix!(1, 3, float, Major.row));
alias Matrix1x3fc = typeof(matrix!(1, 3, float, Major.column));
alias Matrix1x3s = typeof(matrix!(1, 3, short, defaultMajor));
alias Matrix1x3sr = typeof(matrix!(1, 3, short, Major.row));
alias Rvector3s = typeof(matrix!(1, 3, short, Major.row));
alias Matrix1x3sc = typeof(matrix!(1, 3, short, Major.column));
alias Matrix1x3d = typeof(matrix!(1, 3, double, defaultMajor));
alias Matrix1x3dr = typeof(matrix!(1, 3, double, Major.row));
alias Rvector3d = typeof(matrix!(1, 3, double, Major.row));
alias Matrix1x3dc = typeof(matrix!(1, 3, double, Major.column));
alias Matrix1x3cf = typeof(matrix!(1, 3, cfloat, defaultMajor));
alias Matrix1x3cfr = typeof(matrix!(1, 3, cfloat, Major.row));
alias Rvector3cf = typeof(matrix!(1, 3, cfloat, Major.row));
alias Matrix1x3cfc = typeof(matrix!(1, 3, cfloat, Major.column));
alias Matrix1x3ub = typeof(matrix!(1, 3, ubyte, defaultMajor));
alias Matrix1x3ubr = typeof(matrix!(1, 3, ubyte, Major.row));
alias Rvector3ub = typeof(matrix!(1, 3, ubyte, Major.row));
alias Matrix1x3ubc = typeof(matrix!(1, 3, ubyte, Major.column));
alias Matrix1x3b = typeof(matrix!(1, 3, byte, defaultMajor));
alias Matrix1x3br = typeof(matrix!(1, 3, byte, Major.row));
alias Rvector3b = typeof(matrix!(1, 3, byte, Major.row));
alias Matrix1x3bc = typeof(matrix!(1, 3, byte, Major.column));
alias Matrix1x3cr = typeof(matrix!(1, 3, creal, defaultMajor));
alias Matrix1x3crr = typeof(matrix!(1, 3, creal, Major.row));
alias Rvector3cr = typeof(matrix!(1, 3, creal, Major.row));
alias Matrix1x3crc = typeof(matrix!(1, 3, creal, Major.column));
alias Matrix1x3i = typeof(matrix!(1, 3, int, defaultMajor));
alias Matrix1x3ir = typeof(matrix!(1, 3, int, Major.row));
alias Rvector3i = typeof(matrix!(1, 3, int, Major.row));
alias Matrix1x3ic = typeof(matrix!(1, 3, int, Major.column));
alias Matrix1x3ul = typeof(matrix!(1, 3, ulong, defaultMajor));
alias Matrix1x3ulr = typeof(matrix!(1, 3, ulong, Major.row));
alias Rvector3ul = typeof(matrix!(1, 3, ulong, Major.row));
alias Matrix1x3ulc = typeof(matrix!(1, 3, ulong, Major.column));
alias Matrix1x3l = typeof(matrix!(1, 3, long, defaultMajor));
alias Matrix1x3lr = typeof(matrix!(1, 3, long, Major.row));
alias Rvector3l = typeof(matrix!(1, 3, long, Major.row));
alias Matrix1x3lc = typeof(matrix!(1, 3, long, Major.column));
alias Matrix1x3us = typeof(matrix!(1, 3, ushort, defaultMajor));
alias Matrix1x3usr = typeof(matrix!(1, 3, ushort, Major.row));
alias Rvector3us = typeof(matrix!(1, 3, ushort, Major.row));
alias Matrix1x3usc = typeof(matrix!(1, 3, ushort, Major.column));
alias Matrix1x3r = typeof(matrix!(1, 3, real, defaultMajor));
alias Matrix1x3rr = typeof(matrix!(1, 3, real, Major.row));
alias Rvector3r = typeof(matrix!(1, 3, real, Major.row));
alias Matrix1x3rc = typeof(matrix!(1, 3, real, Major.column));
alias Matrix1x3cd = typeof(matrix!(1, 3, cdouble, defaultMajor));
alias Matrix1x3cdr = typeof(matrix!(1, 3, cdouble, Major.row));
alias Rvector3cd = typeof(matrix!(1, 3, cdouble, Major.row));
alias Matrix1x3cdc = typeof(matrix!(1, 3, cdouble, Major.column));
alias Matrix1x4ui = typeof(matrix!(1, 4, uint, defaultMajor));
alias Matrix1x4uir = typeof(matrix!(1, 4, uint, Major.row));
alias Rvector4ui = typeof(matrix!(1, 4, uint, Major.row));
alias Matrix1x4uic = typeof(matrix!(1, 4, uint, Major.column));
alias Matrix1x4f = typeof(matrix!(1, 4, float, defaultMajor));
alias Matrix1x4fr = typeof(matrix!(1, 4, float, Major.row));
alias Rvector4f = typeof(matrix!(1, 4, float, Major.row));
alias Matrix1x4fc = typeof(matrix!(1, 4, float, Major.column));
alias Matrix1x4s = typeof(matrix!(1, 4, short, defaultMajor));
alias Matrix1x4sr = typeof(matrix!(1, 4, short, Major.row));
alias Rvector4s = typeof(matrix!(1, 4, short, Major.row));
alias Matrix1x4sc = typeof(matrix!(1, 4, short, Major.column));
alias Matrix1x4d = typeof(matrix!(1, 4, double, defaultMajor));
alias Matrix1x4dr = typeof(matrix!(1, 4, double, Major.row));
alias Rvector4d = typeof(matrix!(1, 4, double, Major.row));
alias Matrix1x4dc = typeof(matrix!(1, 4, double, Major.column));
alias Matrix1x4cf = typeof(matrix!(1, 4, cfloat, defaultMajor));
alias Matrix1x4cfr = typeof(matrix!(1, 4, cfloat, Major.row));
alias Rvector4cf = typeof(matrix!(1, 4, cfloat, Major.row));
alias Matrix1x4cfc = typeof(matrix!(1, 4, cfloat, Major.column));
alias Matrix1x4ub = typeof(matrix!(1, 4, ubyte, defaultMajor));
alias Matrix1x4ubr = typeof(matrix!(1, 4, ubyte, Major.row));
alias Rvector4ub = typeof(matrix!(1, 4, ubyte, Major.row));
alias Matrix1x4ubc = typeof(matrix!(1, 4, ubyte, Major.column));
alias Matrix1x4b = typeof(matrix!(1, 4, byte, defaultMajor));
alias Matrix1x4br = typeof(matrix!(1, 4, byte, Major.row));
alias Rvector4b = typeof(matrix!(1, 4, byte, Major.row));
alias Matrix1x4bc = typeof(matrix!(1, 4, byte, Major.column));
alias Matrix1x4cr = typeof(matrix!(1, 4, creal, defaultMajor));
alias Matrix1x4crr = typeof(matrix!(1, 4, creal, Major.row));
alias Rvector4cr = typeof(matrix!(1, 4, creal, Major.row));
alias Matrix1x4crc = typeof(matrix!(1, 4, creal, Major.column));
alias Matrix1x4i = typeof(matrix!(1, 4, int, defaultMajor));
alias Matrix1x4ir = typeof(matrix!(1, 4, int, Major.row));
alias Rvector4i = typeof(matrix!(1, 4, int, Major.row));
alias Matrix1x4ic = typeof(matrix!(1, 4, int, Major.column));
alias Matrix1x4ul = typeof(matrix!(1, 4, ulong, defaultMajor));
alias Matrix1x4ulr = typeof(matrix!(1, 4, ulong, Major.row));
alias Rvector4ul = typeof(matrix!(1, 4, ulong, Major.row));
alias Matrix1x4ulc = typeof(matrix!(1, 4, ulong, Major.column));
alias Matrix1x4l = typeof(matrix!(1, 4, long, defaultMajor));
alias Matrix1x4lr = typeof(matrix!(1, 4, long, Major.row));
alias Rvector4l = typeof(matrix!(1, 4, long, Major.row));
alias Matrix1x4lc = typeof(matrix!(1, 4, long, Major.column));
alias Matrix1x4us = typeof(matrix!(1, 4, ushort, defaultMajor));
alias Matrix1x4usr = typeof(matrix!(1, 4, ushort, Major.row));
alias Rvector4us = typeof(matrix!(1, 4, ushort, Major.row));
alias Matrix1x4usc = typeof(matrix!(1, 4, ushort, Major.column));
alias Matrix1x4r = typeof(matrix!(1, 4, real, defaultMajor));
alias Matrix1x4rr = typeof(matrix!(1, 4, real, Major.row));
alias Rvector4r = typeof(matrix!(1, 4, real, Major.row));
alias Matrix1x4rc = typeof(matrix!(1, 4, real, Major.column));
alias Matrix1x4cd = typeof(matrix!(1, 4, cdouble, defaultMajor));
alias Matrix1x4cdr = typeof(matrix!(1, 4, cdouble, Major.row));
alias Rvector4cd = typeof(matrix!(1, 4, cdouble, Major.row));
alias Matrix1x4cdc = typeof(matrix!(1, 4, cdouble, Major.column));
alias Matrix2x1ui = typeof(matrix!(2, 1, uint, defaultMajor));
alias Matrix2x1uir = typeof(matrix!(2, 1, uint, Major.row));
alias Matrix2x1uic = typeof(matrix!(2, 1, uint, Major.column));
alias Cvector2ui = typeof(matrix!(2, 1, uint, Major.column));
alias Matrix2x1f = typeof(matrix!(2, 1, float, defaultMajor));
alias Matrix2x1fr = typeof(matrix!(2, 1, float, Major.row));
alias Matrix2x1fc = typeof(matrix!(2, 1, float, Major.column));
alias Cvector2f = typeof(matrix!(2, 1, float, Major.column));
alias Matrix2x1s = typeof(matrix!(2, 1, short, defaultMajor));
alias Matrix2x1sr = typeof(matrix!(2, 1, short, Major.row));
alias Matrix2x1sc = typeof(matrix!(2, 1, short, Major.column));
alias Cvector2s = typeof(matrix!(2, 1, short, Major.column));
alias Matrix2x1d = typeof(matrix!(2, 1, double, defaultMajor));
alias Matrix2x1dr = typeof(matrix!(2, 1, double, Major.row));
alias Matrix2x1dc = typeof(matrix!(2, 1, double, Major.column));
alias Cvector2d = typeof(matrix!(2, 1, double, Major.column));
alias Matrix2x1cf = typeof(matrix!(2, 1, cfloat, defaultMajor));
alias Matrix2x1cfr = typeof(matrix!(2, 1, cfloat, Major.row));
alias Matrix2x1cfc = typeof(matrix!(2, 1, cfloat, Major.column));
alias Cvector2cf = typeof(matrix!(2, 1, cfloat, Major.column));
alias Matrix2x1ub = typeof(matrix!(2, 1, ubyte, defaultMajor));
alias Matrix2x1ubr = typeof(matrix!(2, 1, ubyte, Major.row));
alias Matrix2x1ubc = typeof(matrix!(2, 1, ubyte, Major.column));
alias Cvector2ub = typeof(matrix!(2, 1, ubyte, Major.column));
alias Matrix2x1b = typeof(matrix!(2, 1, byte, defaultMajor));
alias Matrix2x1br = typeof(matrix!(2, 1, byte, Major.row));
alias Matrix2x1bc = typeof(matrix!(2, 1, byte, Major.column));
alias Cvector2b = typeof(matrix!(2, 1, byte, Major.column));
alias Matrix2x1cr = typeof(matrix!(2, 1, creal, defaultMajor));
alias Matrix2x1crr = typeof(matrix!(2, 1, creal, Major.row));
alias Matrix2x1crc = typeof(matrix!(2, 1, creal, Major.column));
alias Cvector2cr = typeof(matrix!(2, 1, creal, Major.column));
alias Matrix2x1i = typeof(matrix!(2, 1, int, defaultMajor));
alias Matrix2x1ir = typeof(matrix!(2, 1, int, Major.row));
alias Matrix2x1ic = typeof(matrix!(2, 1, int, Major.column));
alias Cvector2i = typeof(matrix!(2, 1, int, Major.column));
alias Matrix2x1ul = typeof(matrix!(2, 1, ulong, defaultMajor));
alias Matrix2x1ulr = typeof(matrix!(2, 1, ulong, Major.row));
alias Matrix2x1ulc = typeof(matrix!(2, 1, ulong, Major.column));
alias Cvector2ul = typeof(matrix!(2, 1, ulong, Major.column));
alias Matrix2x1l = typeof(matrix!(2, 1, long, defaultMajor));
alias Matrix2x1lr = typeof(matrix!(2, 1, long, Major.row));
alias Matrix2x1lc = typeof(matrix!(2, 1, long, Major.column));
alias Cvector2l = typeof(matrix!(2, 1, long, Major.column));
alias Matrix2x1us = typeof(matrix!(2, 1, ushort, defaultMajor));
alias Matrix2x1usr = typeof(matrix!(2, 1, ushort, Major.row));
alias Matrix2x1usc = typeof(matrix!(2, 1, ushort, Major.column));
alias Cvector2us = typeof(matrix!(2, 1, ushort, Major.column));
alias Matrix2x1r = typeof(matrix!(2, 1, real, defaultMajor));
alias Matrix2x1rr = typeof(matrix!(2, 1, real, Major.row));
alias Matrix2x1rc = typeof(matrix!(2, 1, real, Major.column));
alias Cvector2r = typeof(matrix!(2, 1, real, Major.column));
alias Matrix2x1cd = typeof(matrix!(2, 1, cdouble, defaultMajor));
alias Matrix2x1cdr = typeof(matrix!(2, 1, cdouble, Major.row));
alias Matrix2x1cdc = typeof(matrix!(2, 1, cdouble, Major.column));
alias Cvector2cd = typeof(matrix!(2, 1, cdouble, Major.column));
alias Matrix2x2ui = typeof(matrix!(2, 2, uint, defaultMajor));
alias Matrix2ui = typeof(matrix!(2, 2, uint, defaultMajor));
alias Matrix2x2uir = typeof(matrix!(2, 2, uint, Major.row));
alias Matrix2uir = typeof(matrix!(2, 2, uint, Major.row));
alias Matrix2x2uic = typeof(matrix!(2, 2, uint, Major.column));
alias Matrix2uic = typeof(matrix!(2, 2, uint, Major.column));
alias Matrix2x2f = typeof(matrix!(2, 2, float, defaultMajor));
alias Matrix2f = typeof(matrix!(2, 2, float, defaultMajor));
alias Matrix2x2fr = typeof(matrix!(2, 2, float, Major.row));
alias Matrix2fr = typeof(matrix!(2, 2, float, Major.row));
alias Matrix2x2fc = typeof(matrix!(2, 2, float, Major.column));
alias Matrix2fc = typeof(matrix!(2, 2, float, Major.column));
alias Matrix2x2s = typeof(matrix!(2, 2, short, defaultMajor));
alias Matrix2s = typeof(matrix!(2, 2, short, defaultMajor));
alias Matrix2x2sr = typeof(matrix!(2, 2, short, Major.row));
alias Matrix2sr = typeof(matrix!(2, 2, short, Major.row));
alias Matrix2x2sc = typeof(matrix!(2, 2, short, Major.column));
alias Matrix2sc = typeof(matrix!(2, 2, short, Major.column));
alias Matrix2x2d = typeof(matrix!(2, 2, double, defaultMajor));
alias Matrix2d = typeof(matrix!(2, 2, double, defaultMajor));
alias Matrix2x2dr = typeof(matrix!(2, 2, double, Major.row));
alias Matrix2dr = typeof(matrix!(2, 2, double, Major.row));
alias Matrix2x2dc = typeof(matrix!(2, 2, double, Major.column));
alias Matrix2dc = typeof(matrix!(2, 2, double, Major.column));
alias Matrix2x2cf = typeof(matrix!(2, 2, cfloat, defaultMajor));
alias Matrix2cf = typeof(matrix!(2, 2, cfloat, defaultMajor));
alias Matrix2x2cfr = typeof(matrix!(2, 2, cfloat, Major.row));
alias Matrix2cfr = typeof(matrix!(2, 2, cfloat, Major.row));
alias Matrix2x2cfc = typeof(matrix!(2, 2, cfloat, Major.column));
alias Matrix2cfc = typeof(matrix!(2, 2, cfloat, Major.column));
alias Matrix2x2ub = typeof(matrix!(2, 2, ubyte, defaultMajor));
alias Matrix2ub = typeof(matrix!(2, 2, ubyte, defaultMajor));
alias Matrix2x2ubr = typeof(matrix!(2, 2, ubyte, Major.row));
alias Matrix2ubr = typeof(matrix!(2, 2, ubyte, Major.row));
alias Matrix2x2ubc = typeof(matrix!(2, 2, ubyte, Major.column));
alias Matrix2ubc = typeof(matrix!(2, 2, ubyte, Major.column));
alias Matrix2x2b = typeof(matrix!(2, 2, byte, defaultMajor));
alias Matrix2b = typeof(matrix!(2, 2, byte, defaultMajor));
alias Matrix2x2br = typeof(matrix!(2, 2, byte, Major.row));
alias Matrix2br = typeof(matrix!(2, 2, byte, Major.row));
alias Matrix2x2bc = typeof(matrix!(2, 2, byte, Major.column));
alias Matrix2bc = typeof(matrix!(2, 2, byte, Major.column));
alias Matrix2x2cr = typeof(matrix!(2, 2, creal, defaultMajor));
alias Matrix2cr = typeof(matrix!(2, 2, creal, defaultMajor));
alias Matrix2x2crr = typeof(matrix!(2, 2, creal, Major.row));
alias Matrix2crr = typeof(matrix!(2, 2, creal, Major.row));
alias Matrix2x2crc = typeof(matrix!(2, 2, creal, Major.column));
alias Matrix2crc = typeof(matrix!(2, 2, creal, Major.column));
alias Matrix2x2i = typeof(matrix!(2, 2, int, defaultMajor));
alias Matrix2i = typeof(matrix!(2, 2, int, defaultMajor));
alias Matrix2x2ir = typeof(matrix!(2, 2, int, Major.row));
alias Matrix2ir = typeof(matrix!(2, 2, int, Major.row));
alias Matrix2x2ic = typeof(matrix!(2, 2, int, Major.column));
alias Matrix2ic = typeof(matrix!(2, 2, int, Major.column));
alias Matrix2x2ul = typeof(matrix!(2, 2, ulong, defaultMajor));
alias Matrix2ul = typeof(matrix!(2, 2, ulong, defaultMajor));
alias Matrix2x2ulr = typeof(matrix!(2, 2, ulong, Major.row));
alias Matrix2ulr = typeof(matrix!(2, 2, ulong, Major.row));
alias Matrix2x2ulc = typeof(matrix!(2, 2, ulong, Major.column));
alias Matrix2ulc = typeof(matrix!(2, 2, ulong, Major.column));
alias Matrix2x2l = typeof(matrix!(2, 2, long, defaultMajor));
alias Matrix2l = typeof(matrix!(2, 2, long, defaultMajor));
alias Matrix2x2lr = typeof(matrix!(2, 2, long, Major.row));
alias Matrix2lr = typeof(matrix!(2, 2, long, Major.row));
alias Matrix2x2lc = typeof(matrix!(2, 2, long, Major.column));
alias Matrix2lc = typeof(matrix!(2, 2, long, Major.column));
alias Matrix2x2us = typeof(matrix!(2, 2, ushort, defaultMajor));
alias Matrix2us = typeof(matrix!(2, 2, ushort, defaultMajor));
alias Matrix2x2usr = typeof(matrix!(2, 2, ushort, Major.row));
alias Matrix2usr = typeof(matrix!(2, 2, ushort, Major.row));
alias Matrix2x2usc = typeof(matrix!(2, 2, ushort, Major.column));
alias Matrix2usc = typeof(matrix!(2, 2, ushort, Major.column));
alias Matrix2x2r = typeof(matrix!(2, 2, real, defaultMajor));
alias Matrix2r = typeof(matrix!(2, 2, real, defaultMajor));
alias Matrix2x2rr = typeof(matrix!(2, 2, real, Major.row));
alias Matrix2rr = typeof(matrix!(2, 2, real, Major.row));
alias Matrix2x2rc = typeof(matrix!(2, 2, real, Major.column));
alias Matrix2rc = typeof(matrix!(2, 2, real, Major.column));
alias Matrix2x2cd = typeof(matrix!(2, 2, cdouble, defaultMajor));
alias Matrix2cd = typeof(matrix!(2, 2, cdouble, defaultMajor));
alias Matrix2x2cdr = typeof(matrix!(2, 2, cdouble, Major.row));
alias Matrix2cdr = typeof(matrix!(2, 2, cdouble, Major.row));
alias Matrix2x2cdc = typeof(matrix!(2, 2, cdouble, Major.column));
alias Matrix2cdc = typeof(matrix!(2, 2, cdouble, Major.column));
alias Matrix2x3ui = typeof(matrix!(2, 3, uint, defaultMajor));
alias Matrix2x3uir = typeof(matrix!(2, 3, uint, Major.row));
alias Matrix2x3uic = typeof(matrix!(2, 3, uint, Major.column));
alias Matrix2x3f = typeof(matrix!(2, 3, float, defaultMajor));
alias Matrix2x3fr = typeof(matrix!(2, 3, float, Major.row));
alias Matrix2x3fc = typeof(matrix!(2, 3, float, Major.column));
alias Matrix2x3s = typeof(matrix!(2, 3, short, defaultMajor));
alias Matrix2x3sr = typeof(matrix!(2, 3, short, Major.row));
alias Matrix2x3sc = typeof(matrix!(2, 3, short, Major.column));
alias Matrix2x3d = typeof(matrix!(2, 3, double, defaultMajor));
alias Matrix2x3dr = typeof(matrix!(2, 3, double, Major.row));
alias Matrix2x3dc = typeof(matrix!(2, 3, double, Major.column));
alias Matrix2x3cf = typeof(matrix!(2, 3, cfloat, defaultMajor));
alias Matrix2x3cfr = typeof(matrix!(2, 3, cfloat, Major.row));
alias Matrix2x3cfc = typeof(matrix!(2, 3, cfloat, Major.column));
alias Matrix2x3ub = typeof(matrix!(2, 3, ubyte, defaultMajor));
alias Matrix2x3ubr = typeof(matrix!(2, 3, ubyte, Major.row));
alias Matrix2x3ubc = typeof(matrix!(2, 3, ubyte, Major.column));
alias Matrix2x3b = typeof(matrix!(2, 3, byte, defaultMajor));
alias Matrix2x3br = typeof(matrix!(2, 3, byte, Major.row));
alias Matrix2x3bc = typeof(matrix!(2, 3, byte, Major.column));
alias Matrix2x3cr = typeof(matrix!(2, 3, creal, defaultMajor));
alias Matrix2x3crr = typeof(matrix!(2, 3, creal, Major.row));
alias Matrix2x3crc = typeof(matrix!(2, 3, creal, Major.column));
alias Matrix2x3i = typeof(matrix!(2, 3, int, defaultMajor));
alias Matrix2x3ir = typeof(matrix!(2, 3, int, Major.row));
alias Matrix2x3ic = typeof(matrix!(2, 3, int, Major.column));
alias Matrix2x3ul = typeof(matrix!(2, 3, ulong, defaultMajor));
alias Matrix2x3ulr = typeof(matrix!(2, 3, ulong, Major.row));
alias Matrix2x3ulc = typeof(matrix!(2, 3, ulong, Major.column));
alias Matrix2x3l = typeof(matrix!(2, 3, long, defaultMajor));
alias Matrix2x3lr = typeof(matrix!(2, 3, long, Major.row));
alias Matrix2x3lc = typeof(matrix!(2, 3, long, Major.column));
alias Matrix2x3us = typeof(matrix!(2, 3, ushort, defaultMajor));
alias Matrix2x3usr = typeof(matrix!(2, 3, ushort, Major.row));
alias Matrix2x3usc = typeof(matrix!(2, 3, ushort, Major.column));
alias Matrix2x3r = typeof(matrix!(2, 3, real, defaultMajor));
alias Matrix2x3rr = typeof(matrix!(2, 3, real, Major.row));
alias Matrix2x3rc = typeof(matrix!(2, 3, real, Major.column));
alias Matrix2x3cd = typeof(matrix!(2, 3, cdouble, defaultMajor));
alias Matrix2x3cdr = typeof(matrix!(2, 3, cdouble, Major.row));
alias Matrix2x3cdc = typeof(matrix!(2, 3, cdouble, Major.column));
alias Matrix2x4ui = typeof(matrix!(2, 4, uint, defaultMajor));
alias Matrix2x4uir = typeof(matrix!(2, 4, uint, Major.row));
alias Matrix2x4uic = typeof(matrix!(2, 4, uint, Major.column));
alias Matrix2x4f = typeof(matrix!(2, 4, float, defaultMajor));
alias Matrix2x4fr = typeof(matrix!(2, 4, float, Major.row));
alias Matrix2x4fc = typeof(matrix!(2, 4, float, Major.column));
alias Matrix2x4s = typeof(matrix!(2, 4, short, defaultMajor));
alias Matrix2x4sr = typeof(matrix!(2, 4, short, Major.row));
alias Matrix2x4sc = typeof(matrix!(2, 4, short, Major.column));
alias Matrix2x4d = typeof(matrix!(2, 4, double, defaultMajor));
alias Matrix2x4dr = typeof(matrix!(2, 4, double, Major.row));
alias Matrix2x4dc = typeof(matrix!(2, 4, double, Major.column));
alias Matrix2x4cf = typeof(matrix!(2, 4, cfloat, defaultMajor));
alias Matrix2x4cfr = typeof(matrix!(2, 4, cfloat, Major.row));
alias Matrix2x4cfc = typeof(matrix!(2, 4, cfloat, Major.column));
alias Matrix2x4ub = typeof(matrix!(2, 4, ubyte, defaultMajor));
alias Matrix2x4ubr = typeof(matrix!(2, 4, ubyte, Major.row));
alias Matrix2x4ubc = typeof(matrix!(2, 4, ubyte, Major.column));
alias Matrix2x4b = typeof(matrix!(2, 4, byte, defaultMajor));
alias Matrix2x4br = typeof(matrix!(2, 4, byte, Major.row));
alias Matrix2x4bc = typeof(matrix!(2, 4, byte, Major.column));
alias Matrix2x4cr = typeof(matrix!(2, 4, creal, defaultMajor));
alias Matrix2x4crr = typeof(matrix!(2, 4, creal, Major.row));
alias Matrix2x4crc = typeof(matrix!(2, 4, creal, Major.column));
alias Matrix2x4i = typeof(matrix!(2, 4, int, defaultMajor));
alias Matrix2x4ir = typeof(matrix!(2, 4, int, Major.row));
alias Matrix2x4ic = typeof(matrix!(2, 4, int, Major.column));
alias Matrix2x4ul = typeof(matrix!(2, 4, ulong, defaultMajor));
alias Matrix2x4ulr = typeof(matrix!(2, 4, ulong, Major.row));
alias Matrix2x4ulc = typeof(matrix!(2, 4, ulong, Major.column));
alias Matrix2x4l = typeof(matrix!(2, 4, long, defaultMajor));
alias Matrix2x4lr = typeof(matrix!(2, 4, long, Major.row));
alias Matrix2x4lc = typeof(matrix!(2, 4, long, Major.column));
alias Matrix2x4us = typeof(matrix!(2, 4, ushort, defaultMajor));
alias Matrix2x4usr = typeof(matrix!(2, 4, ushort, Major.row));
alias Matrix2x4usc = typeof(matrix!(2, 4, ushort, Major.column));
alias Matrix2x4r = typeof(matrix!(2, 4, real, defaultMajor));
alias Matrix2x4rr = typeof(matrix!(2, 4, real, Major.row));
alias Matrix2x4rc = typeof(matrix!(2, 4, real, Major.column));
alias Matrix2x4cd = typeof(matrix!(2, 4, cdouble, defaultMajor));
alias Matrix2x4cdr = typeof(matrix!(2, 4, cdouble, Major.row));
alias Matrix2x4cdc = typeof(matrix!(2, 4, cdouble, Major.column));
alias Matrix3x1ui = typeof(matrix!(3, 1, uint, defaultMajor));
alias Matrix3x1uir = typeof(matrix!(3, 1, uint, Major.row));
alias Matrix3x1uic = typeof(matrix!(3, 1, uint, Major.column));
alias Cvector3ui = typeof(matrix!(3, 1, uint, Major.column));
alias Matrix3x1f = typeof(matrix!(3, 1, float, defaultMajor));
alias Matrix3x1fr = typeof(matrix!(3, 1, float, Major.row));
alias Matrix3x1fc = typeof(matrix!(3, 1, float, Major.column));
alias Cvector3f = typeof(matrix!(3, 1, float, Major.column));
alias Matrix3x1s = typeof(matrix!(3, 1, short, defaultMajor));
alias Matrix3x1sr = typeof(matrix!(3, 1, short, Major.row));
alias Matrix3x1sc = typeof(matrix!(3, 1, short, Major.column));
alias Cvector3s = typeof(matrix!(3, 1, short, Major.column));
alias Matrix3x1d = typeof(matrix!(3, 1, double, defaultMajor));
alias Matrix3x1dr = typeof(matrix!(3, 1, double, Major.row));
alias Matrix3x1dc = typeof(matrix!(3, 1, double, Major.column));
alias Cvector3d = typeof(matrix!(3, 1, double, Major.column));
alias Matrix3x1cf = typeof(matrix!(3, 1, cfloat, defaultMajor));
alias Matrix3x1cfr = typeof(matrix!(3, 1, cfloat, Major.row));
alias Matrix3x1cfc = typeof(matrix!(3, 1, cfloat, Major.column));
alias Cvector3cf = typeof(matrix!(3, 1, cfloat, Major.column));
alias Matrix3x1ub = typeof(matrix!(3, 1, ubyte, defaultMajor));
alias Matrix3x1ubr = typeof(matrix!(3, 1, ubyte, Major.row));
alias Matrix3x1ubc = typeof(matrix!(3, 1, ubyte, Major.column));
alias Cvector3ub = typeof(matrix!(3, 1, ubyte, Major.column));
alias Matrix3x1b = typeof(matrix!(3, 1, byte, defaultMajor));
alias Matrix3x1br = typeof(matrix!(3, 1, byte, Major.row));
alias Matrix3x1bc = typeof(matrix!(3, 1, byte, Major.column));
alias Cvector3b = typeof(matrix!(3, 1, byte, Major.column));
alias Matrix3x1cr = typeof(matrix!(3, 1, creal, defaultMajor));
alias Matrix3x1crr = typeof(matrix!(3, 1, creal, Major.row));
alias Matrix3x1crc = typeof(matrix!(3, 1, creal, Major.column));
alias Cvector3cr = typeof(matrix!(3, 1, creal, Major.column));
alias Matrix3x1i = typeof(matrix!(3, 1, int, defaultMajor));
alias Matrix3x1ir = typeof(matrix!(3, 1, int, Major.row));
alias Matrix3x1ic = typeof(matrix!(3, 1, int, Major.column));
alias Cvector3i = typeof(matrix!(3, 1, int, Major.column));
alias Matrix3x1ul = typeof(matrix!(3, 1, ulong, defaultMajor));
alias Matrix3x1ulr = typeof(matrix!(3, 1, ulong, Major.row));
alias Matrix3x1ulc = typeof(matrix!(3, 1, ulong, Major.column));
alias Cvector3ul = typeof(matrix!(3, 1, ulong, Major.column));
alias Matrix3x1l = typeof(matrix!(3, 1, long, defaultMajor));
alias Matrix3x1lr = typeof(matrix!(3, 1, long, Major.row));
alias Matrix3x1lc = typeof(matrix!(3, 1, long, Major.column));
alias Cvector3l = typeof(matrix!(3, 1, long, Major.column));
alias Matrix3x1us = typeof(matrix!(3, 1, ushort, defaultMajor));
alias Matrix3x1usr = typeof(matrix!(3, 1, ushort, Major.row));
alias Matrix3x1usc = typeof(matrix!(3, 1, ushort, Major.column));
alias Cvector3us = typeof(matrix!(3, 1, ushort, Major.column));
alias Matrix3x1r = typeof(matrix!(3, 1, real, defaultMajor));
alias Matrix3x1rr = typeof(matrix!(3, 1, real, Major.row));
alias Matrix3x1rc = typeof(matrix!(3, 1, real, Major.column));
alias Cvector3r = typeof(matrix!(3, 1, real, Major.column));
alias Matrix3x1cd = typeof(matrix!(3, 1, cdouble, defaultMajor));
alias Matrix3x1cdr = typeof(matrix!(3, 1, cdouble, Major.row));
alias Matrix3x1cdc = typeof(matrix!(3, 1, cdouble, Major.column));
alias Cvector3cd = typeof(matrix!(3, 1, cdouble, Major.column));
alias Matrix3x2ui = typeof(matrix!(3, 2, uint, defaultMajor));
alias Matrix3x2uir = typeof(matrix!(3, 2, uint, Major.row));
alias Matrix3x2uic = typeof(matrix!(3, 2, uint, Major.column));
alias Matrix3x2f = typeof(matrix!(3, 2, float, defaultMajor));
alias Matrix3x2fr = typeof(matrix!(3, 2, float, Major.row));
alias Matrix3x2fc = typeof(matrix!(3, 2, float, Major.column));
alias Matrix3x2s = typeof(matrix!(3, 2, short, defaultMajor));
alias Matrix3x2sr = typeof(matrix!(3, 2, short, Major.row));
alias Matrix3x2sc = typeof(matrix!(3, 2, short, Major.column));
alias Matrix3x2d = typeof(matrix!(3, 2, double, defaultMajor));
alias Matrix3x2dr = typeof(matrix!(3, 2, double, Major.row));
alias Matrix3x2dc = typeof(matrix!(3, 2, double, Major.column));
alias Matrix3x2cf = typeof(matrix!(3, 2, cfloat, defaultMajor));
alias Matrix3x2cfr = typeof(matrix!(3, 2, cfloat, Major.row));
alias Matrix3x2cfc = typeof(matrix!(3, 2, cfloat, Major.column));
alias Matrix3x2ub = typeof(matrix!(3, 2, ubyte, defaultMajor));
alias Matrix3x2ubr = typeof(matrix!(3, 2, ubyte, Major.row));
alias Matrix3x2ubc = typeof(matrix!(3, 2, ubyte, Major.column));
alias Matrix3x2b = typeof(matrix!(3, 2, byte, defaultMajor));
alias Matrix3x2br = typeof(matrix!(3, 2, byte, Major.row));
alias Matrix3x2bc = typeof(matrix!(3, 2, byte, Major.column));
alias Matrix3x2cr = typeof(matrix!(3, 2, creal, defaultMajor));
alias Matrix3x2crr = typeof(matrix!(3, 2, creal, Major.row));
alias Matrix3x2crc = typeof(matrix!(3, 2, creal, Major.column));
alias Matrix3x2i = typeof(matrix!(3, 2, int, defaultMajor));
alias Matrix3x2ir = typeof(matrix!(3, 2, int, Major.row));
alias Matrix3x2ic = typeof(matrix!(3, 2, int, Major.column));
alias Matrix3x2ul = typeof(matrix!(3, 2, ulong, defaultMajor));
alias Matrix3x2ulr = typeof(matrix!(3, 2, ulong, Major.row));
alias Matrix3x2ulc = typeof(matrix!(3, 2, ulong, Major.column));
alias Matrix3x2l = typeof(matrix!(3, 2, long, defaultMajor));
alias Matrix3x2lr = typeof(matrix!(3, 2, long, Major.row));
alias Matrix3x2lc = typeof(matrix!(3, 2, long, Major.column));
alias Matrix3x2us = typeof(matrix!(3, 2, ushort, defaultMajor));
alias Matrix3x2usr = typeof(matrix!(3, 2, ushort, Major.row));
alias Matrix3x2usc = typeof(matrix!(3, 2, ushort, Major.column));
alias Matrix3x2r = typeof(matrix!(3, 2, real, defaultMajor));
alias Matrix3x2rr = typeof(matrix!(3, 2, real, Major.row));
alias Matrix3x2rc = typeof(matrix!(3, 2, real, Major.column));
alias Matrix3x2cd = typeof(matrix!(3, 2, cdouble, defaultMajor));
alias Matrix3x2cdr = typeof(matrix!(3, 2, cdouble, Major.row));
alias Matrix3x2cdc = typeof(matrix!(3, 2, cdouble, Major.column));
alias Matrix3x3ui = typeof(matrix!(3, 3, uint, defaultMajor));
alias Matrix3ui = typeof(matrix!(3, 3, uint, defaultMajor));
alias Matrix3x3uir = typeof(matrix!(3, 3, uint, Major.row));
alias Matrix3uir = typeof(matrix!(3, 3, uint, Major.row));
alias Matrix3x3uic = typeof(matrix!(3, 3, uint, Major.column));
alias Matrix3uic = typeof(matrix!(3, 3, uint, Major.column));
alias Matrix3x3f = typeof(matrix!(3, 3, float, defaultMajor));
alias Matrix3f = typeof(matrix!(3, 3, float, defaultMajor));
alias Matrix3x3fr = typeof(matrix!(3, 3, float, Major.row));
alias Matrix3fr = typeof(matrix!(3, 3, float, Major.row));
alias Matrix3x3fc = typeof(matrix!(3, 3, float, Major.column));
alias Matrix3fc = typeof(matrix!(3, 3, float, Major.column));
alias Matrix3x3s = typeof(matrix!(3, 3, short, defaultMajor));
alias Matrix3s = typeof(matrix!(3, 3, short, defaultMajor));
alias Matrix3x3sr = typeof(matrix!(3, 3, short, Major.row));
alias Matrix3sr = typeof(matrix!(3, 3, short, Major.row));
alias Matrix3x3sc = typeof(matrix!(3, 3, short, Major.column));
alias Matrix3sc = typeof(matrix!(3, 3, short, Major.column));
alias Matrix3x3d = typeof(matrix!(3, 3, double, defaultMajor));
alias Matrix3d = typeof(matrix!(3, 3, double, defaultMajor));
alias Matrix3x3dr = typeof(matrix!(3, 3, double, Major.row));
alias Matrix3dr = typeof(matrix!(3, 3, double, Major.row));
alias Matrix3x3dc = typeof(matrix!(3, 3, double, Major.column));
alias Matrix3dc = typeof(matrix!(3, 3, double, Major.column));
alias Matrix3x3cf = typeof(matrix!(3, 3, cfloat, defaultMajor));
alias Matrix3cf = typeof(matrix!(3, 3, cfloat, defaultMajor));
alias Matrix3x3cfr = typeof(matrix!(3, 3, cfloat, Major.row));
alias Matrix3cfr = typeof(matrix!(3, 3, cfloat, Major.row));
alias Matrix3x3cfc = typeof(matrix!(3, 3, cfloat, Major.column));
alias Matrix3cfc = typeof(matrix!(3, 3, cfloat, Major.column));
alias Matrix3x3ub = typeof(matrix!(3, 3, ubyte, defaultMajor));
alias Matrix3ub = typeof(matrix!(3, 3, ubyte, defaultMajor));
alias Matrix3x3ubr = typeof(matrix!(3, 3, ubyte, Major.row));
alias Matrix3ubr = typeof(matrix!(3, 3, ubyte, Major.row));
alias Matrix3x3ubc = typeof(matrix!(3, 3, ubyte, Major.column));
alias Matrix3ubc = typeof(matrix!(3, 3, ubyte, Major.column));
alias Matrix3x3b = typeof(matrix!(3, 3, byte, defaultMajor));
alias Matrix3b = typeof(matrix!(3, 3, byte, defaultMajor));
alias Matrix3x3br = typeof(matrix!(3, 3, byte, Major.row));
alias Matrix3br = typeof(matrix!(3, 3, byte, Major.row));
alias Matrix3x3bc = typeof(matrix!(3, 3, byte, Major.column));
alias Matrix3bc = typeof(matrix!(3, 3, byte, Major.column));
alias Matrix3x3cr = typeof(matrix!(3, 3, creal, defaultMajor));
alias Matrix3cr = typeof(matrix!(3, 3, creal, defaultMajor));
alias Matrix3x3crr = typeof(matrix!(3, 3, creal, Major.row));
alias Matrix3crr = typeof(matrix!(3, 3, creal, Major.row));
alias Matrix3x3crc = typeof(matrix!(3, 3, creal, Major.column));
alias Matrix3crc = typeof(matrix!(3, 3, creal, Major.column));
alias Matrix3x3i = typeof(matrix!(3, 3, int, defaultMajor));
alias Matrix3i = typeof(matrix!(3, 3, int, defaultMajor));
alias Matrix3x3ir = typeof(matrix!(3, 3, int, Major.row));
alias Matrix3ir = typeof(matrix!(3, 3, int, Major.row));
alias Matrix3x3ic = typeof(matrix!(3, 3, int, Major.column));
alias Matrix3ic = typeof(matrix!(3, 3, int, Major.column));
alias Matrix3x3ul = typeof(matrix!(3, 3, ulong, defaultMajor));
alias Matrix3ul = typeof(matrix!(3, 3, ulong, defaultMajor));
alias Matrix3x3ulr = typeof(matrix!(3, 3, ulong, Major.row));
alias Matrix3ulr = typeof(matrix!(3, 3, ulong, Major.row));
alias Matrix3x3ulc = typeof(matrix!(3, 3, ulong, Major.column));
alias Matrix3ulc = typeof(matrix!(3, 3, ulong, Major.column));
alias Matrix3x3l = typeof(matrix!(3, 3, long, defaultMajor));
alias Matrix3l = typeof(matrix!(3, 3, long, defaultMajor));
alias Matrix3x3lr = typeof(matrix!(3, 3, long, Major.row));
alias Matrix3lr = typeof(matrix!(3, 3, long, Major.row));
alias Matrix3x3lc = typeof(matrix!(3, 3, long, Major.column));
alias Matrix3lc = typeof(matrix!(3, 3, long, Major.column));
alias Matrix3x3us = typeof(matrix!(3, 3, ushort, defaultMajor));
alias Matrix3us = typeof(matrix!(3, 3, ushort, defaultMajor));
alias Matrix3x3usr = typeof(matrix!(3, 3, ushort, Major.row));
alias Matrix3usr = typeof(matrix!(3, 3, ushort, Major.row));
alias Matrix3x3usc = typeof(matrix!(3, 3, ushort, Major.column));
alias Matrix3usc = typeof(matrix!(3, 3, ushort, Major.column));
alias Matrix3x3r = typeof(matrix!(3, 3, real, defaultMajor));
alias Matrix3r = typeof(matrix!(3, 3, real, defaultMajor));
alias Matrix3x3rr = typeof(matrix!(3, 3, real, Major.row));
alias Matrix3rr = typeof(matrix!(3, 3, real, Major.row));
alias Matrix3x3rc = typeof(matrix!(3, 3, real, Major.column));
alias Matrix3rc = typeof(matrix!(3, 3, real, Major.column));
alias Matrix3x3cd = typeof(matrix!(3, 3, cdouble, defaultMajor));
alias Matrix3cd = typeof(matrix!(3, 3, cdouble, defaultMajor));
alias Matrix3x3cdr = typeof(matrix!(3, 3, cdouble, Major.row));
alias Matrix3cdr = typeof(matrix!(3, 3, cdouble, Major.row));
alias Matrix3x3cdc = typeof(matrix!(3, 3, cdouble, Major.column));
alias Matrix3cdc = typeof(matrix!(3, 3, cdouble, Major.column));
alias Matrix3x4ui = typeof(matrix!(3, 4, uint, defaultMajor));
alias Matrix3x4uir = typeof(matrix!(3, 4, uint, Major.row));
alias Matrix3x4uic = typeof(matrix!(3, 4, uint, Major.column));
alias Matrix3x4f = typeof(matrix!(3, 4, float, defaultMajor));
alias Matrix3x4fr = typeof(matrix!(3, 4, float, Major.row));
alias Matrix3x4fc = typeof(matrix!(3, 4, float, Major.column));
alias Matrix3x4s = typeof(matrix!(3, 4, short, defaultMajor));
alias Matrix3x4sr = typeof(matrix!(3, 4, short, Major.row));
alias Matrix3x4sc = typeof(matrix!(3, 4, short, Major.column));
alias Matrix3x4d = typeof(matrix!(3, 4, double, defaultMajor));
alias Matrix3x4dr = typeof(matrix!(3, 4, double, Major.row));
alias Matrix3x4dc = typeof(matrix!(3, 4, double, Major.column));
alias Matrix3x4cf = typeof(matrix!(3, 4, cfloat, defaultMajor));
alias Matrix3x4cfr = typeof(matrix!(3, 4, cfloat, Major.row));
alias Matrix3x4cfc = typeof(matrix!(3, 4, cfloat, Major.column));
alias Matrix3x4ub = typeof(matrix!(3, 4, ubyte, defaultMajor));
alias Matrix3x4ubr = typeof(matrix!(3, 4, ubyte, Major.row));
alias Matrix3x4ubc = typeof(matrix!(3, 4, ubyte, Major.column));
alias Matrix3x4b = typeof(matrix!(3, 4, byte, defaultMajor));
alias Matrix3x4br = typeof(matrix!(3, 4, byte, Major.row));
alias Matrix3x4bc = typeof(matrix!(3, 4, byte, Major.column));
alias Matrix3x4cr = typeof(matrix!(3, 4, creal, defaultMajor));
alias Matrix3x4crr = typeof(matrix!(3, 4, creal, Major.row));
alias Matrix3x4crc = typeof(matrix!(3, 4, creal, Major.column));
alias Matrix3x4i = typeof(matrix!(3, 4, int, defaultMajor));
alias Matrix3x4ir = typeof(matrix!(3, 4, int, Major.row));
alias Matrix3x4ic = typeof(matrix!(3, 4, int, Major.column));
alias Matrix3x4ul = typeof(matrix!(3, 4, ulong, defaultMajor));
alias Matrix3x4ulr = typeof(matrix!(3, 4, ulong, Major.row));
alias Matrix3x4ulc = typeof(matrix!(3, 4, ulong, Major.column));
alias Matrix3x4l = typeof(matrix!(3, 4, long, defaultMajor));
alias Matrix3x4lr = typeof(matrix!(3, 4, long, Major.row));
alias Matrix3x4lc = typeof(matrix!(3, 4, long, Major.column));
alias Matrix3x4us = typeof(matrix!(3, 4, ushort, defaultMajor));
alias Matrix3x4usr = typeof(matrix!(3, 4, ushort, Major.row));
alias Matrix3x4usc = typeof(matrix!(3, 4, ushort, Major.column));
alias Matrix3x4r = typeof(matrix!(3, 4, real, defaultMajor));
alias Matrix3x4rr = typeof(matrix!(3, 4, real, Major.row));
alias Matrix3x4rc = typeof(matrix!(3, 4, real, Major.column));
alias Matrix3x4cd = typeof(matrix!(3, 4, cdouble, defaultMajor));
alias Matrix3x4cdr = typeof(matrix!(3, 4, cdouble, Major.row));
alias Matrix3x4cdc = typeof(matrix!(3, 4, cdouble, Major.column));
alias Matrix4x1ui = typeof(matrix!(4, 1, uint, defaultMajor));
alias Matrix4x1uir = typeof(matrix!(4, 1, uint, Major.row));
alias Matrix4x1uic = typeof(matrix!(4, 1, uint, Major.column));
alias Cvector4ui = typeof(matrix!(4, 1, uint, Major.column));
alias Matrix4x1f = typeof(matrix!(4, 1, float, defaultMajor));
alias Matrix4x1fr = typeof(matrix!(4, 1, float, Major.row));
alias Matrix4x1fc = typeof(matrix!(4, 1, float, Major.column));
alias Cvector4f = typeof(matrix!(4, 1, float, Major.column));
alias Matrix4x1s = typeof(matrix!(4, 1, short, defaultMajor));
alias Matrix4x1sr = typeof(matrix!(4, 1, short, Major.row));
alias Matrix4x1sc = typeof(matrix!(4, 1, short, Major.column));
alias Cvector4s = typeof(matrix!(4, 1, short, Major.column));
alias Matrix4x1d = typeof(matrix!(4, 1, double, defaultMajor));
alias Matrix4x1dr = typeof(matrix!(4, 1, double, Major.row));
alias Matrix4x1dc = typeof(matrix!(4, 1, double, Major.column));
alias Cvector4d = typeof(matrix!(4, 1, double, Major.column));
alias Matrix4x1cf = typeof(matrix!(4, 1, cfloat, defaultMajor));
alias Matrix4x1cfr = typeof(matrix!(4, 1, cfloat, Major.row));
alias Matrix4x1cfc = typeof(matrix!(4, 1, cfloat, Major.column));
alias Cvector4cf = typeof(matrix!(4, 1, cfloat, Major.column));
alias Matrix4x1ub = typeof(matrix!(4, 1, ubyte, defaultMajor));
alias Matrix4x1ubr = typeof(matrix!(4, 1, ubyte, Major.row));
alias Matrix4x1ubc = typeof(matrix!(4, 1, ubyte, Major.column));
alias Cvector4ub = typeof(matrix!(4, 1, ubyte, Major.column));
alias Matrix4x1b = typeof(matrix!(4, 1, byte, defaultMajor));
alias Matrix4x1br = typeof(matrix!(4, 1, byte, Major.row));
alias Matrix4x1bc = typeof(matrix!(4, 1, byte, Major.column));
alias Cvector4b = typeof(matrix!(4, 1, byte, Major.column));
alias Matrix4x1cr = typeof(matrix!(4, 1, creal, defaultMajor));
alias Matrix4x1crr = typeof(matrix!(4, 1, creal, Major.row));
alias Matrix4x1crc = typeof(matrix!(4, 1, creal, Major.column));
alias Cvector4cr = typeof(matrix!(4, 1, creal, Major.column));
alias Matrix4x1i = typeof(matrix!(4, 1, int, defaultMajor));
alias Matrix4x1ir = typeof(matrix!(4, 1, int, Major.row));
alias Matrix4x1ic = typeof(matrix!(4, 1, int, Major.column));
alias Cvector4i = typeof(matrix!(4, 1, int, Major.column));
alias Matrix4x1ul = typeof(matrix!(4, 1, ulong, defaultMajor));
alias Matrix4x1ulr = typeof(matrix!(4, 1, ulong, Major.row));
alias Matrix4x1ulc = typeof(matrix!(4, 1, ulong, Major.column));
alias Cvector4ul = typeof(matrix!(4, 1, ulong, Major.column));
alias Matrix4x1l = typeof(matrix!(4, 1, long, defaultMajor));
alias Matrix4x1lr = typeof(matrix!(4, 1, long, Major.row));
alias Matrix4x1lc = typeof(matrix!(4, 1, long, Major.column));
alias Cvector4l = typeof(matrix!(4, 1, long, Major.column));
alias Matrix4x1us = typeof(matrix!(4, 1, ushort, defaultMajor));
alias Matrix4x1usr = typeof(matrix!(4, 1, ushort, Major.row));
alias Matrix4x1usc = typeof(matrix!(4, 1, ushort, Major.column));
alias Cvector4us = typeof(matrix!(4, 1, ushort, Major.column));
alias Matrix4x1r = typeof(matrix!(4, 1, real, defaultMajor));
alias Matrix4x1rr = typeof(matrix!(4, 1, real, Major.row));
alias Matrix4x1rc = typeof(matrix!(4, 1, real, Major.column));
alias Cvector4r = typeof(matrix!(4, 1, real, Major.column));
alias Matrix4x1cd = typeof(matrix!(4, 1, cdouble, defaultMajor));
alias Matrix4x1cdr = typeof(matrix!(4, 1, cdouble, Major.row));
alias Matrix4x1cdc = typeof(matrix!(4, 1, cdouble, Major.column));
alias Cvector4cd = typeof(matrix!(4, 1, cdouble, Major.column));
alias Matrix4x2ui = typeof(matrix!(4, 2, uint, defaultMajor));
alias Matrix4x2uir = typeof(matrix!(4, 2, uint, Major.row));
alias Matrix4x2uic = typeof(matrix!(4, 2, uint, Major.column));
alias Matrix4x2f = typeof(matrix!(4, 2, float, defaultMajor));
alias Matrix4x2fr = typeof(matrix!(4, 2, float, Major.row));
alias Matrix4x2fc = typeof(matrix!(4, 2, float, Major.column));
alias Matrix4x2s = typeof(matrix!(4, 2, short, defaultMajor));
alias Matrix4x2sr = typeof(matrix!(4, 2, short, Major.row));
alias Matrix4x2sc = typeof(matrix!(4, 2, short, Major.column));
alias Matrix4x2d = typeof(matrix!(4, 2, double, defaultMajor));
alias Matrix4x2dr = typeof(matrix!(4, 2, double, Major.row));
alias Matrix4x2dc = typeof(matrix!(4, 2, double, Major.column));
alias Matrix4x2cf = typeof(matrix!(4, 2, cfloat, defaultMajor));
alias Matrix4x2cfr = typeof(matrix!(4, 2, cfloat, Major.row));
alias Matrix4x2cfc = typeof(matrix!(4, 2, cfloat, Major.column));
alias Matrix4x2ub = typeof(matrix!(4, 2, ubyte, defaultMajor));
alias Matrix4x2ubr = typeof(matrix!(4, 2, ubyte, Major.row));
alias Matrix4x2ubc = typeof(matrix!(4, 2, ubyte, Major.column));
alias Matrix4x2b = typeof(matrix!(4, 2, byte, defaultMajor));
alias Matrix4x2br = typeof(matrix!(4, 2, byte, Major.row));
alias Matrix4x2bc = typeof(matrix!(4, 2, byte, Major.column));
alias Matrix4x2cr = typeof(matrix!(4, 2, creal, defaultMajor));
alias Matrix4x2crr = typeof(matrix!(4, 2, creal, Major.row));
alias Matrix4x2crc = typeof(matrix!(4, 2, creal, Major.column));
alias Matrix4x2i = typeof(matrix!(4, 2, int, defaultMajor));
alias Matrix4x2ir = typeof(matrix!(4, 2, int, Major.row));
alias Matrix4x2ic = typeof(matrix!(4, 2, int, Major.column));
alias Matrix4x2ul = typeof(matrix!(4, 2, ulong, defaultMajor));
alias Matrix4x2ulr = typeof(matrix!(4, 2, ulong, Major.row));
alias Matrix4x2ulc = typeof(matrix!(4, 2, ulong, Major.column));
alias Matrix4x2l = typeof(matrix!(4, 2, long, defaultMajor));
alias Matrix4x2lr = typeof(matrix!(4, 2, long, Major.row));
alias Matrix4x2lc = typeof(matrix!(4, 2, long, Major.column));
alias Matrix4x2us = typeof(matrix!(4, 2, ushort, defaultMajor));
alias Matrix4x2usr = typeof(matrix!(4, 2, ushort, Major.row));
alias Matrix4x2usc = typeof(matrix!(4, 2, ushort, Major.column));
alias Matrix4x2r = typeof(matrix!(4, 2, real, defaultMajor));
alias Matrix4x2rr = typeof(matrix!(4, 2, real, Major.row));
alias Matrix4x2rc = typeof(matrix!(4, 2, real, Major.column));
alias Matrix4x2cd = typeof(matrix!(4, 2, cdouble, defaultMajor));
alias Matrix4x2cdr = typeof(matrix!(4, 2, cdouble, Major.row));
alias Matrix4x2cdc = typeof(matrix!(4, 2, cdouble, Major.column));
alias Matrix4x3ui = typeof(matrix!(4, 3, uint, defaultMajor));
alias Matrix4x3uir = typeof(matrix!(4, 3, uint, Major.row));
alias Matrix4x3uic = typeof(matrix!(4, 3, uint, Major.column));
alias Matrix4x3f = typeof(matrix!(4, 3, float, defaultMajor));
alias Matrix4x3fr = typeof(matrix!(4, 3, float, Major.row));
alias Matrix4x3fc = typeof(matrix!(4, 3, float, Major.column));
alias Matrix4x3s = typeof(matrix!(4, 3, short, defaultMajor));
alias Matrix4x3sr = typeof(matrix!(4, 3, short, Major.row));
alias Matrix4x3sc = typeof(matrix!(4, 3, short, Major.column));
alias Matrix4x3d = typeof(matrix!(4, 3, double, defaultMajor));
alias Matrix4x3dr = typeof(matrix!(4, 3, double, Major.row));
alias Matrix4x3dc = typeof(matrix!(4, 3, double, Major.column));
alias Matrix4x3cf = typeof(matrix!(4, 3, cfloat, defaultMajor));
alias Matrix4x3cfr = typeof(matrix!(4, 3, cfloat, Major.row));
alias Matrix4x3cfc = typeof(matrix!(4, 3, cfloat, Major.column));
alias Matrix4x3ub = typeof(matrix!(4, 3, ubyte, defaultMajor));
alias Matrix4x3ubr = typeof(matrix!(4, 3, ubyte, Major.row));
alias Matrix4x3ubc = typeof(matrix!(4, 3, ubyte, Major.column));
alias Matrix4x3b = typeof(matrix!(4, 3, byte, defaultMajor));
alias Matrix4x3br = typeof(matrix!(4, 3, byte, Major.row));
alias Matrix4x3bc = typeof(matrix!(4, 3, byte, Major.column));
alias Matrix4x3cr = typeof(matrix!(4, 3, creal, defaultMajor));
alias Matrix4x3crr = typeof(matrix!(4, 3, creal, Major.row));
alias Matrix4x3crc = typeof(matrix!(4, 3, creal, Major.column));
alias Matrix4x3i = typeof(matrix!(4, 3, int, defaultMajor));
alias Matrix4x3ir = typeof(matrix!(4, 3, int, Major.row));
alias Matrix4x3ic = typeof(matrix!(4, 3, int, Major.column));
alias Matrix4x3ul = typeof(matrix!(4, 3, ulong, defaultMajor));
alias Matrix4x3ulr = typeof(matrix!(4, 3, ulong, Major.row));
alias Matrix4x3ulc = typeof(matrix!(4, 3, ulong, Major.column));
alias Matrix4x3l = typeof(matrix!(4, 3, long, defaultMajor));
alias Matrix4x3lr = typeof(matrix!(4, 3, long, Major.row));
alias Matrix4x3lc = typeof(matrix!(4, 3, long, Major.column));
alias Matrix4x3us = typeof(matrix!(4, 3, ushort, defaultMajor));
alias Matrix4x3usr = typeof(matrix!(4, 3, ushort, Major.row));
alias Matrix4x3usc = typeof(matrix!(4, 3, ushort, Major.column));
alias Matrix4x3r = typeof(matrix!(4, 3, real, defaultMajor));
alias Matrix4x3rr = typeof(matrix!(4, 3, real, Major.row));
alias Matrix4x3rc = typeof(matrix!(4, 3, real, Major.column));
alias Matrix4x3cd = typeof(matrix!(4, 3, cdouble, defaultMajor));
alias Matrix4x3cdr = typeof(matrix!(4, 3, cdouble, Major.row));
alias Matrix4x3cdc = typeof(matrix!(4, 3, cdouble, Major.column));
alias Matrix4x4ui = typeof(matrix!(4, 4, uint, defaultMajor));
alias Matrix4ui = typeof(matrix!(4, 4, uint, defaultMajor));
alias Matrix4x4uir = typeof(matrix!(4, 4, uint, Major.row));
alias Matrix4uir = typeof(matrix!(4, 4, uint, Major.row));
alias Matrix4x4uic = typeof(matrix!(4, 4, uint, Major.column));
alias Matrix4uic = typeof(matrix!(4, 4, uint, Major.column));
alias Matrix4x4f = typeof(matrix!(4, 4, float, defaultMajor));
alias Matrix4f = typeof(matrix!(4, 4, float, defaultMajor));
alias Matrix4x4fr = typeof(matrix!(4, 4, float, Major.row));
alias Matrix4fr = typeof(matrix!(4, 4, float, Major.row));
alias Matrix4x4fc = typeof(matrix!(4, 4, float, Major.column));
alias Matrix4fc = typeof(matrix!(4, 4, float, Major.column));
alias Matrix4x4s = typeof(matrix!(4, 4, short, defaultMajor));
alias Matrix4s = typeof(matrix!(4, 4, short, defaultMajor));
alias Matrix4x4sr = typeof(matrix!(4, 4, short, Major.row));
alias Matrix4sr = typeof(matrix!(4, 4, short, Major.row));
alias Matrix4x4sc = typeof(matrix!(4, 4, short, Major.column));
alias Matrix4sc = typeof(matrix!(4, 4, short, Major.column));
alias Matrix4x4d = typeof(matrix!(4, 4, double, defaultMajor));
alias Matrix4d = typeof(matrix!(4, 4, double, defaultMajor));
alias Matrix4x4dr = typeof(matrix!(4, 4, double, Major.row));
alias Matrix4dr = typeof(matrix!(4, 4, double, Major.row));
alias Matrix4x4dc = typeof(matrix!(4, 4, double, Major.column));
alias Matrix4dc = typeof(matrix!(4, 4, double, Major.column));
alias Matrix4x4cf = typeof(matrix!(4, 4, cfloat, defaultMajor));
alias Matrix4cf = typeof(matrix!(4, 4, cfloat, defaultMajor));
alias Matrix4x4cfr = typeof(matrix!(4, 4, cfloat, Major.row));
alias Matrix4cfr = typeof(matrix!(4, 4, cfloat, Major.row));
alias Matrix4x4cfc = typeof(matrix!(4, 4, cfloat, Major.column));
alias Matrix4cfc = typeof(matrix!(4, 4, cfloat, Major.column));
alias Matrix4x4ub = typeof(matrix!(4, 4, ubyte, defaultMajor));
alias Matrix4ub = typeof(matrix!(4, 4, ubyte, defaultMajor));
alias Matrix4x4ubr = typeof(matrix!(4, 4, ubyte, Major.row));
alias Matrix4ubr = typeof(matrix!(4, 4, ubyte, Major.row));
alias Matrix4x4ubc = typeof(matrix!(4, 4, ubyte, Major.column));
alias Matrix4ubc = typeof(matrix!(4, 4, ubyte, Major.column));
alias Matrix4x4b = typeof(matrix!(4, 4, byte, defaultMajor));
alias Matrix4b = typeof(matrix!(4, 4, byte, defaultMajor));
alias Matrix4x4br = typeof(matrix!(4, 4, byte, Major.row));
alias Matrix4br = typeof(matrix!(4, 4, byte, Major.row));
alias Matrix4x4bc = typeof(matrix!(4, 4, byte, Major.column));
alias Matrix4bc = typeof(matrix!(4, 4, byte, Major.column));
alias Matrix4x4cr = typeof(matrix!(4, 4, creal, defaultMajor));
alias Matrix4cr = typeof(matrix!(4, 4, creal, defaultMajor));
alias Matrix4x4crr = typeof(matrix!(4, 4, creal, Major.row));
alias Matrix4crr = typeof(matrix!(4, 4, creal, Major.row));
alias Matrix4x4crc = typeof(matrix!(4, 4, creal, Major.column));
alias Matrix4crc = typeof(matrix!(4, 4, creal, Major.column));
alias Matrix4x4i = typeof(matrix!(4, 4, int, defaultMajor));
alias Matrix4i = typeof(matrix!(4, 4, int, defaultMajor));
alias Matrix4x4ir = typeof(matrix!(4, 4, int, Major.row));
alias Matrix4ir = typeof(matrix!(4, 4, int, Major.row));
alias Matrix4x4ic = typeof(matrix!(4, 4, int, Major.column));
alias Matrix4ic = typeof(matrix!(4, 4, int, Major.column));
alias Matrix4x4ul = typeof(matrix!(4, 4, ulong, defaultMajor));
alias Matrix4ul = typeof(matrix!(4, 4, ulong, defaultMajor));
alias Matrix4x4ulr = typeof(matrix!(4, 4, ulong, Major.row));
alias Matrix4ulr = typeof(matrix!(4, 4, ulong, Major.row));
alias Matrix4x4ulc = typeof(matrix!(4, 4, ulong, Major.column));
alias Matrix4ulc = typeof(matrix!(4, 4, ulong, Major.column));
alias Matrix4x4l = typeof(matrix!(4, 4, long, defaultMajor));
alias Matrix4l = typeof(matrix!(4, 4, long, defaultMajor));
alias Matrix4x4lr = typeof(matrix!(4, 4, long, Major.row));
alias Matrix4lr = typeof(matrix!(4, 4, long, Major.row));
alias Matrix4x4lc = typeof(matrix!(4, 4, long, Major.column));
alias Matrix4lc = typeof(matrix!(4, 4, long, Major.column));
alias Matrix4x4us = typeof(matrix!(4, 4, ushort, defaultMajor));
alias Matrix4us = typeof(matrix!(4, 4, ushort, defaultMajor));
alias Matrix4x4usr = typeof(matrix!(4, 4, ushort, Major.row));
alias Matrix4usr = typeof(matrix!(4, 4, ushort, Major.row));
alias Matrix4x4usc = typeof(matrix!(4, 4, ushort, Major.column));
alias Matrix4usc = typeof(matrix!(4, 4, ushort, Major.column));
alias Matrix4x4r = typeof(matrix!(4, 4, real, defaultMajor));
alias Matrix4r = typeof(matrix!(4, 4, real, defaultMajor));
alias Matrix4x4rr = typeof(matrix!(4, 4, real, Major.row));
alias Matrix4rr = typeof(matrix!(4, 4, real, Major.row));
alias Matrix4x4rc = typeof(matrix!(4, 4, real, Major.column));
alias Matrix4rc = typeof(matrix!(4, 4, real, Major.column));
alias Matrix4x4cd = typeof(matrix!(4, 4, cdouble, defaultMajor));
alias Matrix4cd = typeof(matrix!(4, 4, cdouble, defaultMajor));
alias Matrix4x4cdr = typeof(matrix!(4, 4, cdouble, Major.row));
alias Matrix4cdr = typeof(matrix!(4, 4, cdouble, Major.row));
alias Matrix4x4cdc = typeof(matrix!(4, 4, cdouble, Major.column));
alias Matrix4cdc = typeof(matrix!(4, 4, cdouble, Major.column));