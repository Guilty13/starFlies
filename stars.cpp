#include "stars.h"

#include <vector>
#include <iostream>
#include <cmath>
#include <utility>

#include <stdlib.h>
using namespace std;

template<class T>
ostream& operator <<(ostream& os, const vector<T>& v) {
    for (size_t i = 0; i < v.size(); i++) {
        if (i)
            os << " ";
        os << v[i];
    }
    return os;
}

template<class T>
struct TId {
    using TType = T;
};

template<class TTFloat>
struct TDifferential {
    size_t I;
    TTFloat D;
    // D * dx_I

    TDifferential(size_t i, TTFloat d)
    : I(i)
    , D(d)
    {}
};

template<class TTFloat>
struct TVariation {
    TTFloat X;
    vector<TDifferential<TTFloat>> DX;
    // X + DX; DX is kept sorted by I

    TVariation()
    {}

    TVariation(TTFloat x) // implicit
    : X(x)
    {}

    TVariation(const TDifferential<TTFloat>& dx) // implicit
    : X(TTFloat())
    , DX(1, dx)
    {}

    TVariation(TTFloat x, size_t i, TTFloat d)
    : X(x)
    , DX(1, TDifferential<TTFloat>(i, d))
    {}

    template<class TTOp>
    static vector<TDifferential<TTFloat>> Join(const vector<TDifferential<TTFloat>>& v1, const vector<TDifferential<TTFloat>>& v2, const TTOp& op) {
        vector<TDifferential<TTFloat>> y;
        y.reserve(v1.size() + v2.size());
        size_t j1 = 0, j2 = 0;
        while (j1 < v1.size() || j2 < v2.size()) {
            if (j1 < v1.size() && (j2 >= v2.size() || v1[j1].I < v2[j2].I)) {
                y.push_back(TDifferential<TTFloat>(v1[j1].I, op(v1[j1].D, TTFloat())));
                j1++;
            } else if (j2 < v2.size() && (j1 >= v1.size() || v2[j2].I < v1[j1].I)) {
                y.push_back(TDifferential<TTFloat>(v2[j2].I, op(TTFloat(), v2[j2].D)));
                j2++;
            } else {
                y.push_back(TDifferential<TTFloat>(v1[j1].I, op(v1[j1].D, v2[j2].D)));
                j1++;
                j2++;
            }
        }
        return y;
    }
};

template<class TTFloat>
TDifferential<TTFloat> operator -(const TDifferential<TTFloat>& dx) {
    return {dx.I, -dx.D};
}

template<class TTFloat>
TDifferential<TTFloat> operator *(typename TId<TTFloat>::TType c, const TDifferential<TTFloat>& dx) {
    return {dx.I, c * dx.D};
}

template<class TTFloat>
TDifferential<TTFloat> operator /(const TDifferential<TTFloat>& dx, typename TId<TTFloat>::TType c) {
    return {dx.I, dx.D / c};
}

template<class TTFloat>
TVariation<TTFloat> operator -(const TVariation<TTFloat>& x) {
    TVariation<TTFloat> y(-x.X);
    y.DX.reserve(x.DX.size());
    for (size_t j = 0; j < x.DX.size(); j++) {
        y.DX.push_back(-x.DX[j]);
    }
    return y;
}

template<class TTFloat>
TVariation<TTFloat> operator +(const TVariation<TTFloat>& x1, const TVariation<TTFloat>& x2) {
    TVariation<TTFloat> y(x1.X + x2.X);
    y.DX = TVariation<TTFloat>::Join(x1.DX, x2.DX, [](TTFloat d1, TTFloat d2) {
        return d1 + d2;
    });
    return y;
}

template<class TTFloat>
TVariation<TTFloat> operator +(TTFloat x1, const TVariation<TTFloat>& x2) {
    return TVariation<TTFloat>(x1) + x2;
}


template<class TTFloat>
TVariation<TTFloat>& operator +=(TVariation<TTFloat>& x1, const TVariation<TTFloat>& x2) {
    x1 = x1 + x2;
    return x1;
}

template<class TTFloat>
TVariation<TTFloat> operator -(const TVariation<TTFloat>& x1, const TVariation<TTFloat>& x2) {
    TVariation<TTFloat> y(x1.X - x2.X);
    y.DX = TVariation<TTFloat>::Join(x1.DX, x2.DX, [](TTFloat d1, TTFloat d2) {
        return d1 - d2;
    });
    return y;
}

template<class TTFloat>
TVariation<TTFloat>& operator -=(TVariation<TTFloat>& x1, const TVariation<TTFloat>& x2) {
    x1 = x1 - x2;
    return x1;
}

template<class TTFloat>
TVariation<TTFloat> operator *(const TVariation<TTFloat>& x1, const TVariation<TTFloat>& x2) {
    TVariation<TTFloat> y(x1.X * x2.X);
    y.DX = TVariation<TTFloat>::Join(x1.DX, x2.DX, [&](TTFloat d1, TTFloat d2) {
        return x2.X * d1 + x1.X * d2;
    });
    return y;
}

template<class TTFloat>
TVariation<TTFloat> operator *(TTFloat x1, const TVariation<TTFloat>& x2) {
    return TVariation<TTFloat>(x1) * x2;
}


template<class TTFloat>
TVariation<TTFloat>& operator *=(TVariation<TTFloat>& x1, const TVariation<TTFloat>& x2) {
    x1 = x1 * x2;
    return x1;
}

template<class TTFloat>
TVariation<TTFloat> operator /(const TVariation<TTFloat>& x1, const TVariation<TTFloat>& x2) {
    TTFloat x2Inv = 1. / x2.X;
    TVariation<TTFloat> y(x1.X * x2Inv);
    y.DX = TVariation<TTFloat>::Join(x1.DX, x2.DX, [&](TTFloat d1, TTFloat d2) {
        return d1 * x2Inv - d2 * x1.X * x2Inv * x2Inv;
    });
    return y;
}

template<class TTFloat>
TVariation<TTFloat>& operator /=(TVariation<TTFloat>& x1, const TVariation<TTFloat>& x2) {
    x1 = x1 / x2;
    return x1;
}

template<class TTFloat>
TVariation<TTFloat> exp(const TVariation<TTFloat>& x) {
    TVariation<TTFloat> y(exp(x.X));
    y.DX.reserve(x.DX.size());
    for (size_t j = 0; j < x.DX.size(); j++) {
        y.DX.push_back(y.X * x.DX[j]);
    }
    return y;
}

template<class TTFloat>
TVariation<TTFloat> log(const TVariation<TTFloat>& x) {
    TVariation<TTFloat> y(log(x.X));
    y.DX.reserve(x.DX.size());
    for (size_t i = 0; i < x.DX.size(); i++) {
        y.DX.push_back(x.DX[i] / x.X);
    }
    return y;
}

template<class T>
T pow(T x, T y) {
    return exp(log(x) * y);
}

template<class T>
T sqrt(T x) {
    return pow(x, T(0.5));
}

template<class T>
T logistic(T x) {
    return T(1.) / (T(1.) + exp(-x));
}

template<class T>
T logit(T p) {
    return log(p / (T(1.) - p));
}

struct F {
    template<class T>
    T operator()(const vector<T>& x) {
        T y = T();
        //cout << x.size() << ": " << x << endl;
        for (size_t i = 0; i < x.size() / 2; i++) {
            y += exp(logistic(x[i]) * logistic(x[x.size() - i - 1]));
        }
        //cout << log(y) << endl;
        return log(y);
    }
};

template<class T>
T g1 (const vector<T>& x) {
    return x[0] / x[1] + x[2] / x[3];
}

template<class T>
T g2 (const vector<T>& x) {
    return x[0] * exp(-log(x[1])) + x[2] * exp(-log(x[3]));
}

template<class TTFloat, class TTFunc>
TVariation<TTFloat> NumericDerivative(TTFunc f, const vector<TTFloat>& x0) {
    TVariation<TTFloat> y(f(x0));
    y.DX.reserve(x0.size());
    for (size_t i = 0; i < x0.size(); i++) {
        vector<TTFloat> x = x0;
        x[i] += 1e-3;
        TTFloat y2 = f(x);
        x[i] -= 2e-3;
        TTFloat y1 = f(x);
        y.DX.push_back(TDifferential<TTFloat>(i, (y2-y1) / 2e-3));
    }
    return y;
}

template<class TTFloat>
ostream& operator <<(ostream& os, const TVariation<TTFloat>& x) {
    os << x.X;
    for (size_t j = 0; j < x.DX.size(); j++) {
        os << " + " << x.DX[j].D << "dx_" << x.DX[j].I; 
    }
    return os;
}

template<class TTFloat, class TTFunc>
TVariation<TTFloat> AutoDerivative(TTFunc f, const vector<TTFloat>& x) {
    vector<TVariation<float>> X(x.size());
    for (size_t i = 0; i < x.size(); i++) {
        X[i] = TVariation<TTFloat>(x[i], i, 1);
    }
    return f(X);
}

template <class T>
struct TPointRef {
    // TODO: T x[D];
    //const T& x;
    //const T& y;
    T& x;
    T& y;

    TPointRef<T>& operator =(const TPointRef<T>& p) {
        x = p.x;
        y = p.y;
        return *this;
    }
};

template <class T>
struct TPoint {
    // TODO: T x[D];
    T x;
    T y;

    operator TPointRef<T>() {
        return {x,y};
    }
};

template <class T>
TPoint<T> operator -(TPointRef<T> p) {
    return {-p.x, -p.y};
}

template <class T>
TPoint<T> operator -(TPoint<T> p) {
    return {-p.x, -p.y};
}

template <class T>
TPoint<T> operator *(T c, typename TId<TPointRef<T>>::TType p) {
    return {c * p.x, c * p.y};
}

template <class T>
TPoint<T> operator /(typename TId<TPointRef<T>>::TType p, T c) {
    return {p.x / c, p.y / c};
}

template <class T>
TPoint<T> operator +(TPointRef<T> p1, typename TId<TPointRef<T>>::TType p2) {
    return {p1.x + p2.x, p1.y + p2.y};
}

template <class T>
TPoint<T> operator +(const TPoint<T>& p1, typename TId<TPointRef<T>>::TType p2) {
    return {p1.x + p2.x, p1.y + p2.y};
}

template <class T>
TPoint<T>& operator +=(TPoint<T>& p1, typename TId<TPointRef<T>>::TType p2) {
    p1 = p1 + p2;
    return p1;
}

template <class T>
T operator *(TPointRef<T> p1, typename TId<TPointRef<T>>::TType p2) {
    return p1.x * p2.x + p1.y * p2.y;
}

template <class T>
T operator *(const TPoint<T>& p1, typename TId<TPointRef<T>>::TType p2) {
    return p1.x * p2.x + p1.y * p2.y;
}


template <class T>
struct TLegParamsRef {
    T& ax;
    T& ay;
    T& t;

    TPointRef<T> a = {ax, ay};

    TLegParamsRef<T>& operator =(const TLegParamsRef<T>& l) {
        a = l.a;
        t = l.t;
        return *this;
    }
};

template <class T>
struct TLegParams {
    TPoint<T> a;
    T t;

    operator TLegParamsRef<T>() {
        return {a.x, a.y, t};
    }
};

template <size_t L, class T>
struct TParamsN {
    TLegParamsRef<T> l[L];

    TParamsN(vector<T>& x)
    : TParamsN(x, make_index_sequence<L>())
    {}

    template<size_t... R>
    TParamsN(vector<T>& x, index_sequence<R...>)
    : l({ {x[3*R], x[3*R+1], x[3*R+2]}... })
    {}

    T time() const {
        return _time(make_index_sequence<L>());
    }

    template<size_t... R>
    T _time(index_sequence<R...>) const {
        return (l[R].t + ...);
    }

    static size_t constexpr size() {
        return L;
    }
};

template <class T>
using TParams = TParamsN<7, T>;

struct TG {

    vector<float> X;
    //[px,py,vx,vy]
    float tradeoff = 1.0f;

    TG(const vector<float>& X)
    : X(X)
    {}

    template<class T>
    T error(vector<T>& x) const {
        TParams<T> z {x};
        auto p = P<T>();
        auto v = V<T>();
        for (size_t i = 0; i < z.size(); i++) {
            auto dt = z.l[i].t;
            auto a = z.l[i].a;
            p += dt * v + 0.5f * dt * dt * a;
            v += dt * a;
        }

        return (p*p) + (v*v);
    }

    template<class T>
    TPoint<T> P() const {
        return {X[0], X[1]};
    }

    template<class T>
    TPoint<T> V() const {
        return {X[2], X[3]};
    }

    template<class T>
    T time(vector<T>& _x) const {
        TParams<T> x {_x};
        return x.time();
    }

    template<class T>
    T operator()(vector<T>& x) const {
        
        return - error(x) - tradeoff * time(x);
        /*
        //[ax1,ay1,t1,ax2,ay2,t2];
        auto& px0 = X[0];
        auto& py0 = X[1];
        auto& vx0 = X[2];
        auto& vy0 = X[3];
        auto& ax1 = x[0];
        auto& ay1 = x[1];
        auto& t1 = x[2];
        auto& ax2 = x[3];
        auto& ay2 = x[4];
        auto& t2 = x[5];

        auto px1 = px0 + vx0*t1 + 0.5f*ax1*t1*t1;
        auto py1 = py0 + vy0*t1 + 0.5f*ay1*t1*t1;

        auto vx1 = vx0 + ax1*t1;
        auto vy1 = vy0 + ay1*t1;

        auto px2 = px1 + vx1*t2 + 0.5f*ax2*t2*t2;
        auto py2 = py1 + vy1*t2 + 0.5f*ay2*t2*t2;

        auto vx2 = vx1 + ax2*t2;
        auto vy2 = vy1 + ay2*t2;

        return -(px2*px2 + py2*py2 + vx2*vx2 + vy2*vy2) - 0.1f*(t1 + t2);
        */
    }
};


void FixConstraints(vector<float>& y) { // params...
    TParams<float> p{y};

    for (size_t i = 0; i < p.size(); i++) {
        auto& l = p.l[i];
        auto d2 = l.ax * l.ax + l.ay * l.ay;
        if (d2 > 1.) {
            auto d = sqrt(d2);
            l.ax /= d;
            l.ay /= d;
        }
        if (l.t < 0.) {
            l.t = 0.f;
        }
    }
}

vector<float> Init(const TG& G) {
    size_t legs = TParams<float>::size();
    vector<float> x(legs * 3);

    auto p = G.P<float>();
    auto v = G.V<float>();

    vector<TLegParams<float>> l(3);

    // stage 1: reduce initial speed

    float vL = sqrt(v * v);
    auto a1 = l[0].a = -v / vL;
    auto dt1 = l[0].t = vL;
    
    p += dt1 * v + 0.5f * dt1 * dt1 * a1;
    v += dt1 * a1; // 0.0

    // stage 2: accelerate half way to destination point

    float pL = sqrt(p * p);
    auto a2 = l[1].a = -p / pL;
    auto dt2 = l[1].t = sqrt(pL);

    // stage 3: brake other half way to destination point

    l[2].a = -a2;
    l[2].t = dt2;

    // distribute stages over legs
    TParams<float> z{x};

    float dt = (dt2 * 2.0 + dt1) / z.size();

    for (size_t i = 0, j = 0; ; ) {
        if (z.size() - j == l.size() - i) {
            for (; i < l.size(); i++, j++) {
                z.l[j] = l[i];
            }
            break;
        }

        if (l[i].t < dt) {
            z.l[j] = l[i];
            i++;
            j++;
        } else {
            z.l[j].a = l[i].a;
            z.l[j].t = dt;
            l[i].t -= dt;
            j++;
        }
    }

    return x;
}

vector<float> Compute(const vector<float>& params, bool verbose) {
    TG G(params);

    vector<float> x = Init(G);

    for (size_t m = 0; m < 30; m++) {
        for (size_t n = 0; n < 100; n++) {
            auto g = AutoDerivative(G, x);
            //cout << g << endl;

            float s = 1e-3;
            vector<float> y;
            vector<float> x0 = x;
            float r = G(x);

            while (true) {
                y = x0;
                for (auto& d: g.DX) {
                    y[d.I] += s * d.D;
                }

                FixConstraints(y);

                if (r < G(y)) {
                    r = G(y);
                    x = y;
                    s *= 2;
                }
                else {
                    break;
                }
                if (verbose) {
                    cout << x << ": " << G.time(x) << ": " << G.error(x) << endl;
                }
            }
        }
        G.tradeoff /= 2.0f;
    }
    return x;
}

int main2() {
    //TG G({-50,-40,2,-3});
    float x = getenv("X") ? atof(getenv("X")) : 0.0f;
    float y = getenv("Y") ? atof(getenv("Y")) : 0.0f;
    float vx = getenv("VX") ? atof(getenv("VX")) : -1.0f;
    float vy = getenv("VY") ? atof(getenv("VY")) : 1.0f;

    Compute({x, y, vx, vy}, true);
    return 0;
}

int main1() {
    vector<float> x(6);
    for (size_t i = 0; i < x.size(); i++) {
        x[i] = (i+1);
    }
    vector<TVariation<float>> X(x.size());
    for (size_t i = 0; i < x.size(); i++) {
        X[i] = TVariation<float>(x[i], i, 0.5);
    }

    cout << AutoDerivative(F(), x) << endl << NumericDerivative(F(), x) << endl;
    //cout << f(X) << endl << NumericDerivative(f<float>, x) << endl;
    return 0;
}
