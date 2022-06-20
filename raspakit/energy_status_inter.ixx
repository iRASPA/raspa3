export module energy_status_inter;

import <string>;
import <iostream>;
import <sstream>;
import <cmath>;

export struct EnergyInter
{
    double VanDerWaals;
    double VanDerWaalsTailCorrection;
    double CoulombicReal;
    double CoulombicFourier;
    double totalInter;

    EnergyInter() : VanDerWaals(0.0), VanDerWaalsTailCorrection(0.0), CoulombicReal(0.0), CoulombicFourier(0.0), totalInter(0.0)
    {
    }

    void zero()
    {
        VanDerWaals = 0.0;
        VanDerWaalsTailCorrection = 0.0;
        CoulombicReal = 0.0;
        CoulombicFourier = 0.0;
        totalInter = 0.0;
    }

    void sumTotal()
    {
      totalInter = VanDerWaals + VanDerWaalsTailCorrection + CoulombicReal + CoulombicFourier;
    }

    inline double total() const
    {
        return VanDerWaals + VanDerWaalsTailCorrection + CoulombicReal + CoulombicFourier;
    }

    inline EnergyInter& operator+=(const EnergyInter& b)
    {
        this->VanDerWaals += b.VanDerWaals;
        this->VanDerWaalsTailCorrection += b.VanDerWaalsTailCorrection;
        this->CoulombicReal += b.CoulombicReal;
        this->CoulombicFourier += b.CoulombicFourier;
        this->totalInter += b.totalInter;
        return *this;
    }

    inline EnergyInter& operator-=(const EnergyInter& b)
    {
        this->VanDerWaals -= b.VanDerWaals;
        this->VanDerWaalsTailCorrection -= b.VanDerWaalsTailCorrection;
        this->CoulombicReal -= b.CoulombicReal;
        this->CoulombicFourier -= b.CoulombicFourier;
        this->totalInter -= b.totalInter;
        return *this;
    }

    inline EnergyInter operator-() const
    {
        EnergyInter v;
        v.VanDerWaals = -VanDerWaals;
        v.VanDerWaalsTailCorrection = -VanDerWaalsTailCorrection;
        v.CoulombicReal = -CoulombicReal;
        v.CoulombicFourier = -CoulombicFourier;
        v.totalInter = -totalInter;
        return v;
    }
};

export inline EnergyInter operator+(const EnergyInter& a, const EnergyInter& b)
{
    EnergyInter m{};
    m.VanDerWaals = a.VanDerWaals + b.VanDerWaals;
    m.VanDerWaalsTailCorrection = a.VanDerWaalsTailCorrection + b.VanDerWaalsTailCorrection;
    m.CoulombicReal = a.CoulombicReal + b.CoulombicReal;
    m.CoulombicFourier = a.CoulombicFourier + b.CoulombicFourier;
    m.totalInter = a.totalInter + b.totalInter;
    return m;
}

export inline EnergyInter operator-(const EnergyInter& a, const EnergyInter& b)
{
    EnergyInter m{};
    m.VanDerWaals = a.VanDerWaals - b.VanDerWaals;
    m.VanDerWaalsTailCorrection = a.VanDerWaalsTailCorrection - b.VanDerWaalsTailCorrection;
    m.CoulombicReal = a.CoulombicReal - b.CoulombicReal;
    m.CoulombicFourier = a.CoulombicFourier - b.CoulombicFourier;
    m.totalInter = a.totalInter - b.totalInter;
    return m;
}

export inline EnergyInter operator*(const EnergyInter& a, const EnergyInter& b)
{
    EnergyInter m{};
    m.VanDerWaals = a.VanDerWaals * b.VanDerWaals;
    m.VanDerWaalsTailCorrection = a.VanDerWaalsTailCorrection * b.VanDerWaalsTailCorrection;
    m.CoulombicReal = a.CoulombicReal * b.CoulombicReal;
    m.CoulombicFourier = a.CoulombicFourier * b.CoulombicFourier;
    m.totalInter = a.totalInter * b.totalInter;
    return m;
}

export inline EnergyInter operator*(const double& a, const EnergyInter& b)
{
    EnergyInter m{};
    m.VanDerWaals = a * b.VanDerWaals;
    m.VanDerWaalsTailCorrection = a * b.VanDerWaalsTailCorrection;
    m.CoulombicReal = a * b.CoulombicReal;
    m.CoulombicFourier = a * b.CoulombicFourier;
    m.totalInter = a * b.totalInter;
    return m;
}

export inline EnergyInter operator/(const EnergyInter& a, const double& b)
{
    EnergyInter m{};
    m.VanDerWaals = a.VanDerWaals / b;
    m.VanDerWaalsTailCorrection = a.VanDerWaalsTailCorrection / b;
    m.CoulombicReal = a.CoulombicReal / b;
    m.CoulombicFourier = a.CoulombicFourier / b;
    m.totalInter = a.totalInter / b;
    return m;
}

export inline EnergyInter sqrt(const EnergyInter & a)
{
    EnergyInter m{};
    m.VanDerWaals = std::sqrt(a.VanDerWaals);
    m.VanDerWaalsTailCorrection = std::sqrt(a.VanDerWaalsTailCorrection);
    m.CoulombicReal = std::sqrt(a.CoulombicReal);
    m.CoulombicFourier = std::sqrt(a.CoulombicFourier);
    m.totalInter = std::sqrt(a.totalInter);
    return m;
}

