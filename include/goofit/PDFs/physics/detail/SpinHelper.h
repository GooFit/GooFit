/*
04/22/2016 Christoph Hasse
GPU Adaptation of classes provided in the MINT2 package by Jonas Rademacker.
They are needed for the calculation of the spin factors.
DISCLAIMER:
This code is not sufficiently tested yet and still under heavy development!
*/

#pragma once

#include <goofit/PDFs/physics/DalitzPlotHelpers.h>

namespace GooFit {

class __align__(16) gpuLVec {
  private:
    fptype X;
    fptype Y;
    fptype Z;
    fptype E;

  public:
    __device__ gpuLVec() { X = Y = Z = E = 0; };
    __device__ gpuLVec(fptype x, fptype y, fptype z, fptype e)
        : X(x)
        , Y(y)
        , Z(z)
        , E(e){};

    __device__ void SetXYZE(fptype x, fptype y, fptype z, fptype e) {
        X = x;
        Y = y;
        Z = z;
        E = e;
    }

    __device__ auto Dot(const gpuLVec &rhs) const->fptype { return E * rhs.E - X * rhs.X - Y * rhs.Y - Z * rhs.Z; }
    __device__ inline auto GetX() const->fptype { return X; }
    __device__ inline auto GetY() const->fptype { return Y; }
    __device__ inline auto GetZ() const->fptype { return Z; }
    __device__ inline auto GetE() const->fptype { return E; }

    __device__ inline auto Mag2() const->fptype { return this->Dot(*this); }
    __device__ inline auto M() const->fptype { return sqrt(this->Mag2()); }
    __device__ inline void SetX(fptype a) { X = a; }
    __device__ inline void SetY(fptype a) { Y = a; }
    __device__ inline void SetZ(fptype a) { Z = a; }
    __device__ inline void SetE(fptype a) { E = a; }

    __device__ auto operator+=(const gpuLVec &rhs)->gpuLVec & {
        X += rhs.X;
        Y += rhs.Y;
        Z += rhs.Z;
        E += rhs.E;
        return *this;
    }
    __device__ auto operator-=(const gpuLVec &rhs)->gpuLVec & {
        X -= rhs.X;
        Y -= rhs.Y;
        Z -= rhs.Z;
        E -= rhs.E;
        return *this;
    }
    __device__ auto operator*=(const fptype rhs)->gpuLVec & {
        X *= rhs;
        Y *= rhs;
        Z *= rhs;
        E *= rhs;
        return *this;
    }
};
__device__ auto operator+(gpuLVec lhs, const gpuLVec &rhs) -> gpuLVec { return lhs += rhs; }
__device__ auto operator-(gpuLVec lhs, const gpuLVec &rhs) -> gpuLVec { return lhs -= rhs; }
__device__ auto operator*(gpuLVec lhs, fptype rhs) -> gpuLVec { return lhs *= rhs; }
__device__ auto operator*(fptype lhs, gpuLVec rhs) -> gpuLVec { return rhs *= lhs; }

class ZTspin1 : public gpuLVec {
  public:
    // in decay D -> AB:   q = p(A) - p(B), p = p(A) + p(B)
    __device__ ZTspin1(const gpuLVec &q, const gpuLVec &p, fptype mR)
        : gpuLVec(q - q.Dot(p) * p * (1. / (mR * mR))) {}

    __device__ auto Contract(const gpuLVec &rhs) const -> fptype { return this->Dot(rhs); }
};

class SpinSumV { // spin sum for Vector->PP
  protected:
    gpuLVec _p; // pA + pB (pA, pB are the 4-mom of dgtrs)
    fptype _mR; // nominal mass of resonance.
  public:
    __device__ SpinSumV(const gpuLVec &p, fptype mR)
        : _p(p)
        , _mR(mR) {}

    __device__ auto Dot(const gpuLVec &rhs) const -> gpuLVec { return -1.0 * rhs + _p * (_p.Dot(rhs) / (_mR * _mR)); }
    __device__ auto Sandwich(const gpuLVec &lhs, const gpuLVec &rhs) const -> fptype { return lhs.Dot(this->Dot(rhs)); }
};

class SpinSumT {
  protected:
    SpinSumV _sv;

  public:
    __device__ SpinSumT(const gpuLVec &p, fptype mR)
        : _sv(p, mR) {}
    __device__ auto Sandwich(const gpuLVec &lm, const gpuLVec &ln, const gpuLVec &ra, const gpuLVec &rb) -> fptype {
        fptype manb = _sv.Sandwich(lm, ra) * _sv.Sandwich(ln, rb);
        fptype mbna = _sv.Sandwich(lm, rb) * _sv.Sandwich(ln, ra);
        fptype mnab = _sv.Sandwich(lm, ln) * _sv.Sandwich(ra, rb);

        return (1. / 2.) * (manb + mbna) - (1. / 3.) * mnab;
    }
};

class LorentzMatrix {
  protected:
    gpuLVec _v[4];
    // we'll follow the x, y, z, E convention, i.e. E is 4
    __device__ auto makeZero() -> bool {
        X().SetXYZE(0, 0, 0, 0);
        Y().SetXYZE(0, 0, 0, 0);
        Z().SetXYZE(0, 0, 0, 0);
        E().SetXYZE(0, 0, 0, 0);
        return true;
    }

  public:
    __device__ auto v(int i) const -> const gpuLVec & { return _v[i]; }

    __device__ LorentzMatrix() = default;

    __device__ LorentzMatrix(const gpuLVec p[4]) {
        for(int i = 0; i < 4; i++)
            _v[i] = p[i];
    }
    __device__ LorentzMatrix(const LorentzMatrix &other) {
        for(int i = 0; i < 4; i++)
            _v[i] = other._v[i];
    }
    __device__ auto X() const -> const gpuLVec & { return _v[0]; }
    __device__ auto Y() const -> const gpuLVec & { return _v[1]; }
    __device__ auto Z() const -> const gpuLVec & { return _v[2]; }
    __device__ auto E() const -> const gpuLVec & { return _v[3]; }

    __device__ auto X() -> gpuLVec & { return _v[0]; }
    __device__ auto Y() -> gpuLVec & { return _v[1]; }
    __device__ auto Z() -> gpuLVec & { return _v[2]; }
    __device__ auto E() -> gpuLVec & { return _v[3]; }

    __device__ auto operator[](int i) const -> const gpuLVec & { return _v[i]; }
    __device__ auto operator[](int i) -> gpuLVec & { return _v[i]; }

    __device__ auto add(const LorentzMatrix &other) -> LorentzMatrix & {
        for(int i = 0; i < 4; i++)
            _v[i] += other._v[i];

        return *this;
    }
    __device__ auto subtract(const LorentzMatrix &other) -> LorentzMatrix & {
        for(int i = 0; i < 4; i++)
            _v[i] -= other._v[i];

        return *this;
    }
    __device__ auto mult(fptype s) -> LorentzMatrix & {
        for(auto &i : _v)
            i *= s;

        return *this;
    }
    __device__ auto div(fptype s) -> LorentzMatrix & {
        for(auto &i : _v)
            i *= (1. / s);

        return *this;
    }

    __device__ auto operator+=(const LorentzMatrix &rhs) -> LorentzMatrix & { return add(rhs); }
    __device__ auto operator*=(fptype rhs) -> LorentzMatrix & { return mult(rhs); }
    __device__ auto operator-=(const LorentzMatrix &rhs) -> LorentzMatrix & { return subtract(rhs); }
    __device__ auto operator/=(fptype rhs) -> LorentzMatrix & { return div(rhs); }
    __device__ auto operator=(const LorentzMatrix &other) -> LorentzMatrix & {
        for(int i = 0; i < 4; i++)
            _v[i] = other._v[i];

        return *this;
    }
    __device__ auto operator+(const LorentzMatrix &rhs) const -> LorentzMatrix {
        LorentzMatrix returnVal(*this);
        returnVal += rhs;
        return returnVal;
    }
    __device__ auto operator-(const LorentzMatrix &rhs) const -> LorentzMatrix {
        LorentzMatrix returnVal(*this);
        returnVal -= rhs;
        return returnVal;
    }
    __device__ auto operator*(fptype rhs) const -> LorentzMatrix {
        LorentzMatrix returnVal(*this);
        returnVal *= rhs;
        return returnVal;
    }
    __device__ auto operator/(fptype rhs) const -> LorentzMatrix {
        LorentzMatrix returnVal(*this);
        returnVal /= rhs;
        return returnVal;
    }
};

class SymmLorentzMatrix : public LorentzMatrix {
  protected:
    // SymmLorentzMatrix __gmunu;

    // we'll follow the x, y, z, E convention, i.e. E is 4
    __device__ auto symmetrize() -> bool {
        // clumsy but save
        X().SetY(Y().GetX());
        X().SetZ(Z().GetX());
        X().SetE(E().GetX());

        Y().SetX(X().GetY());
        Y().SetZ(Z().GetY());
        Y().SetE(E().GetY());

        Z().SetX(X().GetZ());
        Z().SetY(Y().GetZ());
        Z().SetE(E().GetZ());

        E().SetX(X().GetE());
        E().SetY(Y().GetE());
        E().SetZ(Z().GetE());
        return true;
    }
    __device__ auto makeZero() -> bool {
        X().SetXYZE(0, 0, 0, 0);
        Y().SetXYZE(0, 0, 0, 0);
        Z().SetXYZE(0, 0, 0, 0);
        E().SetXYZE(0, 0, 0, 0);
        return true;
    }

  public:
    __device__ void makeGmunu() {
        X().SetXYZE(-1, 0, 0, 0);
        Y().SetXYZE(0, -1, 0, 0);
        Z().SetXYZE(0, 0, -1, 0);
        E().SetXYZE(0, 0, 0, 1);
    }
    __device__ auto gmunu() -> const SymmLorentzMatrix &;
    __device__ SymmLorentzMatrix()
        : LorentzMatrix() {}
    __device__ SymmLorentzMatrix(const gpuLVec p[4])
        : LorentzMatrix(p){};

    __device__ SymmLorentzMatrix(const gpuLVec p) {
        X().SetX(p.GetX() * p.GetX());
        X().SetY(p.GetX() * p.GetY());
        X().SetZ(p.GetX() * p.GetZ());
        X().SetE(p.GetX() * p.GetE());

        Y().SetX(p.GetY() * p.GetX());
        Y().SetY(p.GetY() * p.GetY());
        Y().SetZ(p.GetY() * p.GetZ());
        Y().SetE(p.GetY() * p.GetE());

        Z().SetX(p.GetZ() * p.GetX());
        Z().SetY(p.GetZ() * p.GetY());
        Z().SetZ(p.GetZ() * p.GetZ());
        Z().SetE(p.GetZ() * p.GetE());

        E().SetX(p.GetE() * p.GetX());
        E().SetY(p.GetE() * p.GetY());
        E().SetZ(p.GetE() * p.GetZ());
        E().SetE(p.GetE() * p.GetE());
    }
    __device__ SymmLorentzMatrix(const SymmLorentzMatrix &other) = default;

    __device__ auto add(const SymmLorentzMatrix &other) -> SymmLorentzMatrix & {
        for(int i = 0; i < 4; i++)
            _v[i] += other._v[i];

        return *this;
    }
    __device__ auto subtract(const SymmLorentzMatrix &other) -> SymmLorentzMatrix & {
        for(int i = 0; i < 4; i++)
            _v[i] -= other._v[i];

        return *this;
    }
    __device__ auto mult(fptype s) -> SymmLorentzMatrix & {
        for(auto &i : _v)
            i *= s;

        return *this;
    }
    __device__ auto div(fptype s) -> SymmLorentzMatrix & {
        for(auto &i : _v)
            i *= (1. / s);

        return *this;
    }
    __device__ auto Contract(const gpuLVec &vec) -> gpuLVec {
        // M^{mu nu} g_{nu alpha} v^{alpha}
        // M^{mu nu} v_{alpha}
        return vec.GetE() * E() - vec.GetX() * X() - vec.GetY() * Y() - vec.GetZ() * Z();
    }
    __device__ auto Contract_1(const SymmLorentzMatrix &M) -> LorentzMatrix {
        // One pair of indices gets contracted. Since
        // both matrices are symmetric, it doesn't matter which.
        //
        // O^{mu alpha} g_{alpha beta} M^{beta nu} = R^{mu nu}
        // O^{mu alpha} M_{beta}^{nu}
        LorentzMatrix R;
        R.X() = this->Contract(M.X());
        R.Y() = this->Contract(M.Y());
        R.Z() = this->Contract(M.Z());
        R.E() = this->Contract(M.E()); // signs?

        return R;
    }
    __device__ auto Contract_2(const SymmLorentzMatrix &M) -> fptype {
        // Both pairs of indices are contracted.
        // since the matrices are symmetric, it does
        // not matter which index from this with which form M.
        //
        // O^{mu alpha} g_{alpha beta} M^{beta nu} g_{mu nu}
        LorentzMatrix R(Contract_1(M));
        // R^{mu nu} R_{mu nu}
        fptype xx = R.X().GetX(); // Dot(R.X());
        fptype yy = R.Y().GetY(); // Dot(R.Y());
        fptype zz = R.Z().GetZ(); // Dot(R.Z());
        fptype tt = R.E().GetE(); // Dot(R.E());

        return tt - xx - yy - zz; // signs?
    }

    __device__ auto operator+=(const SymmLorentzMatrix &rhs) -> SymmLorentzMatrix & { return add(rhs); }
    __device__ auto operator*=(fptype rhs) -> SymmLorentzMatrix & { return mult(rhs); }
    __device__ auto operator-=(const SymmLorentzMatrix &rhs) -> SymmLorentzMatrix & { return subtract(rhs); }
    __device__ auto operator/=(fptype rhs) -> SymmLorentzMatrix & { return div(rhs); }
    __device__ auto operator=(const SymmLorentzMatrix &other) -> SymmLorentzMatrix & {
        for(int i = 0; i < 4; i++)
            _v[i] = other._v[i];

        return *this;
    }
    __device__ auto operator+(const SymmLorentzMatrix &rhs) const -> SymmLorentzMatrix {
        SymmLorentzMatrix returnVal(*this);
        returnVal += rhs;
        return returnVal;
    }
    __device__ auto operator-(const SymmLorentzMatrix &rhs) const -> SymmLorentzMatrix {
        SymmLorentzMatrix returnVal(*this);
        returnVal -= rhs;
        return returnVal;
    }
    __device__ auto operator*(fptype rhs) const -> SymmLorentzMatrix {
        SymmLorentzMatrix returnVal(*this);
        returnVal *= rhs;
        return returnVal;
    }
    __device__ auto operator/(fptype rhs) const -> SymmLorentzMatrix {
        SymmLorentzMatrix returnVal(*this);
        returnVal /= rhs;
        return returnVal;
    }
};

// __device__ SymmLorentzMatrix operator*(fptype lhs, const SymmLorentzMatrix& rhs);
// __device__ SymmLorentzMatrix operator/(fptype lhs, const SymmLorentzMatrix& rhs);

__device__ auto operator*(fptype lhs, const SymmLorentzMatrix &rhs) -> SymmLorentzMatrix {
    SymmLorentzMatrix returnVal(rhs);
    returnVal *= lhs;
    return returnVal;
}
__device__ auto operator/(fptype lhs, const SymmLorentzMatrix &rhs) -> SymmLorentzMatrix {
    SymmLorentzMatrix returnVal(rhs);
    returnVal /= lhs;
    return returnVal;
}

class ZTspin2 : public SymmLorentzMatrix {
  public:
    // in decay D -> AB:   q = p(A) - p(B), p = p(A) + p(B)
    __device__ ZTspin2(const gpuLVec &q, const gpuLVec &p, fptype mR) {
        ZTspin1 t(q, p, mR);
        SymmLorentzMatrix tt(t);
        SymmLorentzMatrix uu(p);
        uu /= (mR * mR);
        auto gmunu = SymmLorentzMatrix();
        gmunu.makeGmunu();
        // printf("%f %f %f %f\n",gmunu.X().GetX(), gmunu.Y().GetY(), gmunu.Z().GetZ(), gmunu.E().GetE() );
        *this = tt - (1. / 3.) * t.Mag2() * (gmunu - uu);
        // eq 6 in PhysRevD.51.2247
    }

    __device__ auto operator=(const SymmLorentzMatrix &other) -> ZTspin2 & {
        for(int i = 0; i < 4; i++)
            _v[i] = other.v(i);

        return *this;
    }
};

__device__ auto LeviCivita(const gpuLVec &p1, const gpuLVec &p2, const gpuLVec &p3, const gpuLVec &p4) -> fptype {
    // this calculates the determinant of the 4x4 matrix build out of p1,p2,p3,p4
    return p1.GetZ() * p2.GetY() * p3.GetX() * p4.GetE() - p1.GetY() * p2.GetZ() * p3.GetX() * p4.GetE()
           - p1.GetZ() * p2.GetX() * p3.GetY() * p4.GetE() + p1.GetX() * p2.GetZ() * p3.GetY() * p4.GetE()
           + p1.GetY() * p2.GetX() * p3.GetZ() * p4.GetE() - p1.GetX() * p2.GetY() * p3.GetZ() * p4.GetE()
           - p1.GetZ() * p2.GetY() * p3.GetE() * p4.GetX() + p1.GetY() * p2.GetZ() * p3.GetE() * p4.GetX()
           + p1.GetZ() * p2.GetE() * p3.GetY() * p4.GetX() - p1.GetE() * p2.GetZ() * p3.GetY() * p4.GetX()
           - p1.GetY() * p2.GetE() * p3.GetZ() * p4.GetX() + p1.GetE() * p2.GetY() * p3.GetZ() * p4.GetX()
           + p1.GetZ() * p2.GetX() * p3.GetE() * p4.GetY() - p1.GetX() * p2.GetZ() * p3.GetE() * p4.GetY()
           - p1.GetZ() * p2.GetE() * p3.GetX() * p4.GetY() + p1.GetE() * p2.GetZ() * p3.GetX() * p4.GetY()
           + p1.GetX() * p2.GetE() * p3.GetZ() * p4.GetY() - p1.GetE() * p2.GetX() * p3.GetZ() * p4.GetY()
           - p1.GetY() * p2.GetX() * p3.GetE() * p4.GetZ() + p1.GetX() * p2.GetY() * p3.GetE() * p4.GetZ()
           + p1.GetY() * p2.GetE() * p3.GetX() * p4.GetZ() - p1.GetE() * p2.GetY() * p3.GetX() * p4.GetZ()
           - p1.GetX() * p2.GetE() * p3.GetY() * p4.GetZ() + p1.GetE() * p2.GetX() * p3.GetY() * p4.GetZ();
}

__device__ inline auto LeviCivita(const gpuLVec &a, const gpuLVec &b, const gpuLVec &c) -> gpuLVec {
    gpuLVec v;

    v.SetE(-1. * a.GetX() * (b.GetY() * c.GetZ() - b.GetZ() * c.GetY())
           - a.GetY() * (b.GetZ() * c.GetX() - b.GetX() * c.GetZ())
           - a.GetZ() * (b.GetX() * c.GetY() - b.GetY() * c.GetX()));

    v.SetZ(a.GetX() * (b.GetE() * c.GetY() - b.GetY() * c.GetE())
           + a.GetY() * (b.GetX() * c.GetE() - b.GetE() * c.GetX())
           + a.GetE() * (b.GetY() * c.GetX() - b.GetX() * c.GetY()));

    v.SetY(a.GetX() * (b.GetZ() * c.GetE() - b.GetE() * c.GetZ())
           + a.GetZ() * (b.GetE() * c.GetX() - b.GetX() * c.GetE())
           + a.GetE() * (b.GetX() * c.GetZ() - b.GetZ() * c.GetX()));

    v.SetX(a.GetY() * (b.GetE() * c.GetZ() - b.GetZ() * c.GetE())
           + a.GetZ() * (b.GetY() * c.GetE() - b.GetE() * c.GetY())
           + a.GetE() * (b.GetZ() * c.GetY() - b.GetY() * c.GetZ()));

    return v;
}

} // namespace GooFit
