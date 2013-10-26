module gigue.quaternion;

import gigue.core;

import std.format,
       std.functional,
       std.math,
       std.traits;

debug = Quaternion;


/**
四元数
*/
Quaternion!(CommonType!(A, B, C, D)) quaternion(A, B, C, D)(A a, B b, C c, D d) pure nothrow @safe
if(is(Quaternion!(CommonType!(A, B, C, D))))
{
    typeof(return) dst;
    dst.a = a;
    dst.b = b;
    dst.c = c;
    dst.d = d;

    return dst;
}


///
Quaternion!A quaternion(A)(A a) pure nothrow @safe
if(is(Quaternion!A))
{
    typeof(return) dst;
    dst.a = a;

    return dst;
}


/// ditto
Quaternion!(CommonType!(R, ElementType!V)) quaternion(R, V)(R r, V v)
if(is(Quaternion!(CommonType!(R, ElementType!V))))
{
    typeof(return) dst;
    dst.s = r;
    dst.v = v;

    return dst;
}


/// ditto
Quaternion!E quaternion(E)(E[] arr)
if(is(Quaternion!E))
in{
    assert(arr.length == 4);
}
body{
    typeof(return) dst;
    dst._vec4.array[] = arr[];

    return dst;
}



/// ditto
struct Quaternion(S)
if(isScalar!S)
{
    this(E)(Quaternion!E q)
    if(is(E : S))
    {
        this = q;
    }


    this()(SCVector!(int, 4) m)
    {
        this._vec4 = m;
    }


    /// 
    ref inout(S) opIndex(size_t i) pure nothrow @safe inout
    in{
        assert(i < 4);
    }
    body{
        return _vec4[i];
    }


    ref inout(S) s() pure nothrow @safe @property inout { return _vec4[0]; }
    ref inout(S) i() pure nothrow @safe @property inout { return _vec4[1]; }
    ref inout(S) j() pure nothrow @safe @property inout { return _vec4[2]; }
    ref inout(S) k() pure nothrow @safe @property inout { return _vec4[3]; }


    alias a = s;
    alias b = i;
    alias c = j;
    alias d = k;


    @property
    auto v()() pure nothrow @safe inout
    {
        return _vec4.reference.swizzle.bcd;
    }


    void v(V)(V v) @property
    {
        foreach(i; 0 .. 3)
            this._vec4[i+1] = v[i];
    }


    Quaternion!(CommonType!(S, E)) opBinary(string op : "+", E)(in Quaternion!E q) const
    if(!is(CommonType!(S, E) == void))
    {
        return typeof(return)(typeof(typeof(return).init._vec4)(this._vec4 + q._vec4));
    }


    Quaternion!(CommonType!(S, E)) opBinary(string op : "-", E)(in Quaternion!E q) const
    if(!is(CommonType!(S, E) == void))
    {
        return typeof(return)(typeof(typeof(return).init._vec4)(this._vec4 - q._vec4));
    }


    Quaternion!(CommonType!(S, E)) opBinary(string op : "*", E)(in Quaternion!E q) const
    if(!is(CommonType!(S, E) == void))
    {
        return quaternion(this.s * q.s - this.v.dot(q.v), (this.s * q.v) + (q.s * this.v) + (this.v.cross(q.v)));
    }


    auto opBinary(string op : "/", E)(in Quaternion!E q) const
    if(!is(CommonType!(S, E) == void))
    {
        return (this * q.conj) / q.abs ^^ 2;
    }


    Quaternion!(CommonType!(S, E)) opBinary(string op : "+", E)(in E s) const
    if(!is(CommonType!(S, E) == void))
    {
        typeof(return) dst;
        dst = this;
        dst.a += s;
        return dst;
    }


    Quaternion!(CommonType!(S, E)) opBinary(string op  : "-", E)(in E s) const
    if(!is(CommonType!(S, E) == void))
    {
        typeof(return) dst;
        dst = this;
        dst.a -= s;
        return dst;
    }


    Quaternion!(CommonType!(S, E)) opBinary(string op : "*", E)(in E s) const
    if(!is(CommonType!(S, E) == void))
    {
        typeof(return) dst;
        dst = this;
        dst._vec4 *= s;
        return dst;
    }


    Quaternion!(CommonType!(S, E)) opBinary(string op : "/", E)(in E s) const
    if(!is(CommonType!(S, E) == void))
    {
        typeof(return) dst;
        dst = this;
        dst._vec4 /= s;
        return dst;
    }


    Quaternion!(CommonType!(S, E)) opBinaryRight(string op : "+", E)(in E s) const
    if(!is(CommonType!(S, E) == void))
    {
        typeof(return) dst;
        dst = this;
        dst.a += s;
        return dst;
    }


    Quaternion!(CommonType!(S, E)) opBinaryRight(string op : "-", E)(in E s) const
    if(!is(CommonType!(S, E) == void))
    {
        return quaternion!(CommonType!(S, E))(s) - this;
    }


    Quaternion!(CommonType!(S, E)) opBinaryRight(string op : "*", E)(in E s) const
    if(!is(CommonType!(S, E) == void))
    {
        typeof(return) dst;
        dst = this;
        dst._vec4 *= s;
        return dst;
    }


    auto opBinaryRight(string op : "/", E)(in E s) const
    if(!is(CommonType!(S, E) == void))
    {
        return s / this.abs ^^ 2 * this.conj;
    }


    void opAssign(E)(in Quaternion!E q)
    if(is(E : S))
    {
        this._vec4 = q._vec4;
    }


    void opAssign(E)(in E s)
    if(is(E : S))
    {
        this._vec4 = 0;
        this.a = s;
    }


    void opOpAssign(string op, E)(Quaternion!E q)
    if(!is(CommonType!(S, E) == void))
    {
        this = mixin("this " ~ op ~ " q");
    }


    void opOpAssign(string op, E)(E s)
    if(is(E : S))
    {
        this = mixin("this " ~ op ~ " s");
    }


    void toString(scope void delegate(const(char)[]) sink, string formatString) const
    {
        formattedWrite(sink, formatString, _vec4.array);
    }


    bool opEquals(E)(auto ref const Quaternion!E q) pure nothrow @safe const
    {
        foreach(i; 0 .. 4)
            if(this[i] != q[i])
                return false;
        return true;
    }


  private:
    SCVector!(S, 4) _vec4;
}


/// 
unittest {
    // 1 = [1; (0, 0, 0)]な四元数の作成
    auto q = quaternion(1);

    // 添字によるアクセス
    assert(q[0] == 1);
    assert(q[1] == 0);
    assert(q[2] == 0);
    assert(q[3] == 0);


    // 1 + 2i + 3j + 4k = [1; (2, 3, 4)]な四元数の作成
    q = quaternion(1, 2, 3, 4);
    assert(q[0] == 1);
    assert(q[1] == 2);
    assert(q[2] == 3);
    assert(q[3] == 4);

    // a, b, c, dによるアクセス
    assert(q.a == 1);
    assert(q.b == 2);
    assert(q.c == 3);
    assert(q.d == 4);

    // スカラー部であるs, ベクトル部であるvによるアクセス
    assert(q.s == 1);
    assert(q.v == [2, 3, 4].matrix!(3, 1));

    // v = (i, j, k)
    assert(q.i == 2);
    assert(q.j == 3);
    assert(q.k == 4);

    // opIndexやa, b, c, d, i, j, k, s, vへは代入可能
    q.s = 7;
    assert(q[0] == 7);

    // vはベクトルなので、ベクトルを代入可能
    q.v = [4, 5, 6].matrix!(3, 1);
    assert(q[1] == 4);
    assert(q[2] == 5);
    assert(q[3] == 6);

    // スカラー部とベクトル部による四元数の作成
    q = quaternion(8, [9, 10, 11].matrix!(3, 1));
    assert(q[0] == 8);
    assert(q[1] == 9);
    assert(q[2] == 10);
    assert(q[3] == 11);


    // 和
    q = quaternion(1, 2, 3, 4) + quaternion(2, 2, 2, 2);
    assert(q == quaternion(3, 4, 5, 6));

    q = q + 3;
    assert(q == quaternion(6, 4, 5, 6));

    q = 3 + q;
    assert(q == quaternion(9, 4, 5, 6));

    // 複合代入和
    q += q;
    assert(q == quaternion(18, 8, 10, 12));

    q += 3;
    assert(q == quaternion(21, 8, 10, 12));


    // 差
    q = quaternion(1, 2, 3, 4) - quaternion(2, 2, 2, 2);
    assert(q == quaternion(-1, 0, 1, 2));

    q = q - 3;
    assert(q == quaternion(-4, 0, 1, 2));

    q = 3 - q;
    assert(q == quaternion(7, 0, -1, -2));

    // 複合代入和
    q -= q;
    assert(q == quaternion(0, 0, 0, 0));

    q -= 3;
    assert(q == quaternion(-3, 0, 0, 0));


    // 積
    q = quaternion(1, 2, 3, 4) * quaternion(7, 6, 7, 8);
    assert(q == quaternion(-58, 16, 36, 32));

    q = quaternion(1, 2, 3, 4) * 4;
    assert(q == quaternion(4, 8, 12, 16));

    q = 4 * quaternion(1, 2, 3, 4);
    assert(q == quaternion(4, 8, 12, 16));

    q = quaternion(1, 2, 3, 4);
    q *= quaternion(7, 6, 7, 8);
    assert(q == quaternion(-58, 16, 36, 32));

    q = quaternion(1, 2, 3, 4);
    q *= 4;
    assert(q == quaternion(4, 8, 12, 16));


    // 商
    assert((quaternion(-58, 16, 36, 32) / quaternion(7, 6, 7, 8)).approxEqual(quaternion(1, 2, 3, 4)));
    assert(quaternion(4, 8, 12, 16) / 4 == quaternion(1, 2, 3, 4));
    assert((16 / quaternion(1, 2, 3, 4)).approxEqual(quaternion(16) / quaternion(1, 2, 3, 4)));
    auto p = quaternion(-58.0, 16, 36, 32);
    p /= quaternion(7, 6, 7, 8);
    assert(p.approxEqual(quaternion(1, 2, 3, 4)));

    p = quaternion(4, 8, 12, 16);
    p /= 4;
    assert(p.approxEqual(quaternion(1, 2, 3, 4)));
}


/**
四元数の絶対値を返します
*/
auto abs(E)(Quaternion!E q)
{
  static if(isFloatingPoint!E)
    return sqrt(q.a ^^ 2 + q.b ^^ 2 + q.c ^^ 2 + q.d ^^ 2);
  else
    return sqrt(cast(real)q.a ^^ 2 + q.b ^^ 2 + q.c ^^ 2 + q.d ^^ 2);
}


/**
四元数の共役を返します
*/
Quaternion!E conj(E)(const Quaternion!E q) pure nothrow @safe
{
    typeof(return) dst;
    dst.s = q.s;
    dst.v = q.v * -1;
    return dst;
}


/**
approxEqualの四元数バージョン
*/
bool approxEqual(alias pred = std.math.approxEqual, Q1, Q2)(Q1 q1, Q2 q2)
{
    foreach(i; 0 .. 4)
        if(!binaryFun!pred(q1[i], q2[i]))
            return false;
    return true;
}