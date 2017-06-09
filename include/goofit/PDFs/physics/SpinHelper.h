/*
04/22/2016 Christoph Hasse
GPU Adaptation of classes provided in the MINT2 package by Jonas Rademacker.
They are needed for the calculation of the spin factors.
DISCLAIMER:
This code is not sufficently tested yet and still under heavy development!
*/

#pragma once

#include "goofit/PDFs/physics/DalitzPlotHelpers.h"

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

    __device__ fptype Dot(const gpuLVec &rhs) const { return E * rhs.E - X * rhs.X - Y * rhs.Y - Z * rhs.Z; }
    __device__ inline fptype GetX() const { return X; }
    __device__ inline fptype GetY() const { return Y; }
    __device__ inline fptype GetZ() const { return Z; }
    __device__ inline fptype GetE() const { return E; }

    __device__ inline fptype Mag2() const { return this->Dot(*this); }
    __device__ inline fptype M() const { return sqrt(this->Mag2()); }
    __device__ inline void SetX(fptype a) { X = a; }
    __device__ inline void SetY(fptype a) { Y = a; }
    __device__ inline void SetZ(fptype a) { Z = a; }
    __device__ inline void SetE(fptype a) { E = a; }

    __device__ gpuLVec &operator+=(const gpuLVec &rhs) {
        X += rhs.X;
        Y += rhs.Y;
        Z += rhs.Z;
        E += rhs.E;
        return *this;
    }
    __device__ gpuLVec &operator-=(const gpuLVec &rhs) {
        X -= rhs.X;
        Y -= rhs.Y;
        Z -= rhs.Z;
        E -= rhs.E;
        return *this;
    }
    __device__ gpuLVec &operator*=(const fptype rhs) {
        X *= rhs;
        Y *= rhs;
        Z *= rhs;
        E *= rhs;
        return *this;
    }
};
__device__ gpuLVec operator+(gpuLVec lhs, const gpuLVec &rhs) { return lhs += rhs; }
__device__ gpuLVec operator-(gpuLVec lhs, const gpuLVec &rhs) { return lhs -= rhs; }
__device__ gpuLVec operator*(gpuLVec lhs, fptype rhs) { return lhs *= rhs; }
__device__ gpuLVec operator*(fptype lhs, gpuLVec rhs) { return rhs *= lhs; }

class ZTspin1 : public gpuLVec {
  public:
    // in decay D -> AB:   q = p(A) - p(B), p = p(A) + p(B)
    __device__ ZTspin1(const gpuLVec &q, const gpuLVec &p, fptype mR)
        : gpuLVec(q - q.Dot(p) * p * (1. / (mR * mR))) {}

    __device__ fptype Contract(const gpuLVec &rhs) const { return this->Dot(rhs); }
};

class SpinSumV { // spin sum for Vector->PP
  protected:
    gpuLVec _p; // pA + pB (pA, pB are the 4-mom of dgtrs)
    fptype _mR; // nominal mass of resonance.
  public:
    __device__ SpinSumV(const gpuLVec &p, fptype mR)
        : _p(p)
        , _mR(mR) {}

    __device__ gpuLVec Dot(const gpuLVec &rhs) const { return -1.0 * rhs + _p * (_p.Dot(rhs) / (_mR * _mR)); }
    __device__ fptype Sandwich(const gpuLVec &lhs, const gpuLVec &rhs) const { return lhs.Dot(this->Dot(rhs)); }
};

class SpinSumT {
  protected:
    SpinSumV _sv;

  public:
    __device__ SpinSumT(const gpuLVec &p, fptype mR)
        : _sv(p, mR) {}
    __device__ fptype Sandwich(const gpuLVec &lm, const gpuLVec &ln, const gpuLVec &ra, const gpuLVec &rb) {
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
    __device__ bool makeZero() {
        X().SetXYZE(0, 0, 0, 0);
        Y().SetXYZE(0, 0, 0, 0);
        Z().SetXYZE(0, 0, 0, 0);
        E().SetXYZE(0, 0, 0, 0);
        return true;
    }

  public:
    __device__ const gpuLVec &v(int i) const { return _v[i]; }

    __device__ LorentzMatrix() = default;

    __device__ LorentzMatrix(const gpuLVec p[4]) {
        for(int i = 0; i < 4; i++)
            _v[i] = p[i];
    }
    __device__ LorentzMatrix(const LorentzMatrix &other) {
        for(int i = 0; i < 4; i++)
            _v[i] = other._v[i];
    }
    __device__ const gpuLVec &X() const { return _v[0]; }
    __device__ const gpuLVec &Y() const { return _v[1]; }
    __device__ const gpuLVec &Z() const { return _v[2]; }
    __device__ const gpuLVec &E() const { return _v[3]; }

    __device__ gpuLVec &X() { return _v[0]; }
    __device__ gpuLVec &Y() { return _v[1]; }
    __device__ gpuLVec &Z() { return _v[2]; }
    __device__ gpuLVec &E() { return _v[3]; }

    __device__ const gpuLVec &operator[](int i) const { return _v[i]; }
    __device__ gpuLVec &operator[](int i) { return _v[i]; }

    __device__ LorentzMatrix &add(const LorentzMatrix &other) {
        for(int i = 0; i < 4; i++)
            _v[i] += other._v[i];

        return *this;
    }
    __device__ LorentzMatrix &subtract(const LorentzMatrix &other) {
        for(int i = 0; i < 4; i++)
            _v[i] -= other._v[i];

        return *this;
    }
    __device__ LorentzMatrix &mult(fptype s) {
        for(auto &i : _v)
            i *= s;

        return *this;
    }
    __device__ LorentzMatrix &div(fptype s) {
        for(auto &i : _v)
            i *= (1. / s);

        return *this;
    }

    __device__ LorentzMatrix &operator+=(const LorentzMatrix &rhs) { return add(rhs); }
    __device__ LorentzMatrix &operator*=(fptype rhs) { return mult(rhs); }
    __device__ LorentzMatrix &operator-=(const LorentzMatrix &rhs) { return subtract(rhs); }
    __device__ LorentzMatrix &operator/=(fptype rhs) { return div(rhs); }
    __device__ LorentzMatrix &operator=(const LorentzMatrix &other) {
        for(int i = 0; i < 4; i++)
            _v[i] = other._v[i];

        return *this;
    }
    __device__ LorentzMatrix operator+(const LorentzMatrix &rhs) const {
        LorentzMatrix returnVal(*this);
        returnVal += rhs;
        return returnVal;
    }
    __device__ LorentzMatrix operator-(const LorentzMatrix &rhs) const {
        LorentzMatrix returnVal(*this);
        returnVal -= rhs;
        return returnVal;
    }
    __device__ LorentzMatrix operator*(fptype rhs) const {
        LorentzMatrix returnVal(*this);
        returnVal *= rhs;
        return returnVal;
    }
    __device__ LorentzMatrix operator/(fptype rhs) const {
        LorentzMatrix returnVal(*this);
        returnVal /= rhs;
        return returnVal;
    }
};

class SymmLorentzMatrix : public LorentzMatrix {
  protected:
    // SymmLorentzMatrix __gmunu;

    // we'll follow the x, y, z, E convention, i.e. E is 4
    __device__ bool symmetrize() {
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
    __device__ bool makeZero() {
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
    __device__ const SymmLorentzMatrix &gmunu();
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

    __device__ SymmLorentzMatrix &add(const SymmLorentzMatrix &other) {
        for(int i = 0; i < 4; i++)
            _v[i] += other._v[i];

        return *this;
    }
    __device__ SymmLorentzMatrix &subtract(const SymmLorentzMatrix &other) {
        for(int i = 0; i < 4; i++)
            _v[i] -= other._v[i];

        return *this;
    }
    __device__ SymmLorentzMatrix &mult(fptype s) {
        for(auto &i : _v)
            i *= s;

        return *this;
    }
    __device__ SymmLorentzMatrix &div(fptype s) {
        for(auto &i : _v)
            i *= (1. / s);

        return *this;
    }
    __device__ gpuLVec Contract(const gpuLVec &vec) {
        // M^{mu nu} g_{nu alpha} v^{alpha}
        // M^{mu nu} v_{alpha}
        return vec.GetE() * E() - vec.GetX() * X() - vec.GetY() * Y() - vec.GetZ() * Z();
    }
    __device__ LorentzMatrix Contract_1(const SymmLorentzMatrix &M) {
        // One pair of indices gets contracted. Since
        // both matrices are symmetric, it doesnt matter which.
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
    __device__ fptype Contract_2(const SymmLorentzMatrix &M) {
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

    __device__ SymmLorentzMatrix &operator+=(const SymmLorentzMatrix &rhs) { return add(rhs); }
    __device__ SymmLorentzMatrix &operator*=(fptype rhs) { return mult(rhs); }
    __device__ SymmLorentzMatrix &operator-=(const SymmLorentzMatrix &rhs) { return subtract(rhs); }
    __device__ SymmLorentzMatrix &operator/=(fptype rhs) { return div(rhs); }
    __device__ SymmLorentzMatrix &operator=(const SymmLorentzMatrix &other) {
        for(int i = 0; i < 4; i++)
            _v[i] = other._v[i];

        return *this;
    }
    __device__ SymmLorentzMatrix operator+(const SymmLorentzMatrix &rhs) const {
        SymmLorentzMatrix returnVal(*this);
        returnVal += rhs;
        return returnVal;
    }
    __device__ SymmLorentzMatrix operator-(const SymmLorentzMatrix &rhs) const {
        SymmLorentzMatrix returnVal(*this);
        returnVal -= rhs;
        return returnVal;
    }
    __device__ SymmLorentzMatrix operator*(fptype rhs) const {
        SymmLorentzMatrix returnVal(*this);
        returnVal *= rhs;
        return returnVal;
    }
    __device__ SymmLorentzMatrix operator/(fptype rhs) const {
        SymmLorentzMatrix returnVal(*this);
        returnVal /= rhs;
        return returnVal;
    }
};

// __device__ SymmLorentzMatrix operator*(fptype lhs, const SymmLorentzMatrix& rhs);
// __device__ SymmLorentzMatrix operator/(fptype lhs, const SymmLorentzMatrix& rhs);

__device__ SymmLorentzMatrix operator*(fptype lhs, const SymmLorentzMatrix &rhs) {
    SymmLorentzMatrix returnVal(rhs);
    returnVal *= lhs;
    return returnVal;
}
__device__ SymmLorentzMatrix operator/(fptype lhs, const SymmLorentzMatrix &rhs) {
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

    __device__ ZTspin2 &operator=(const SymmLorentzMatrix &other) {
        for(int i = 0; i < 4; i++)
            _v[i] = other.v(i);

        return *this;
    }
};

__device__ fptype LeviCivita(const gpuLVec &p1, const gpuLVec &p2, const gpuLVec &p3, const gpuLVec &p4) {
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

__device__ inline gpuLVec LeviCivita(const gpuLVec &a, const gpuLVec &b, const gpuLVec &c) {
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
