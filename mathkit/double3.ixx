export module double3;

import <cmath>;
import <ostream>;

#if defined(WIN32)
    import <intrin.h>;
#else
    import <immintrin.h>;
#endif


import int3;
import bool3;
import hashcombine;


export union double3
{
    __m256d vec;
    double v[4];
    struct { double x, y, z, w; };

    double3(__m256d const& v): vec(v) {}
    double3& operator = (__m256d const& x) { vec = x; return *this; }
//#ifdef __AVX__
//    double3& load(double const* p) { vec = _mm256_load_pd(p); return *this; }
//#endif
    double3(double x = 0, double y = 0, double z = 0): x(x), y(y), z(z), w(0.0) {}
    double3(int3 a) :x(double(a.x)), y(double(a.y)), z(double(a.z)), w(0.0) {}

    inline double& operator [] (size_t i) { return v[i]; }
    inline const double& operator [] (size_t i) const { return v[i]; }

    double3 normalise();
    double3 fract();
    static double3 randomVectorOnUnitSphere();

    inline double length() { return sqrt(x * x + y * y + z * z); }
    inline double length_squared() const { return (x * x + y * y + z * z); }
    inline double length() const { return sqrt(x * x + y * y + z * z); }
    inline static double3 abs(double3 v1) { return double3(std::abs(v1.x), std::abs(v1.y), std::abs(v1.z)); }
    inline static double3 rint(double3 const& v) { return double3(std::rint(v.x), std::rint(v.y), std::rint(v.z)); }
    inline static double3 fract(double3 const& v) { return double3(v.x - floor(v.x), v.y - floor(v.y), v.z - floor(v.z)); }
    inline static double dot(const double3& v1, const double3& v2) 
    { 
//#ifdef __AVX__
//        __m256d dp = _mm256_mul_pd(v1.vec, v2.vec);
//        __m128d a = _mm256_extractf128_pd(dp, 0);
//        __m128d b = _mm256_extractf128_pd(dp, 1);
//        __m128d c = _mm_add_pd(a, b);
//        __m128d yy = _mm_unpackhi_pd(c, c);
//        __m128d dotproduct = _mm_add_pd(c, yy);
//        return _mm_cvtsd_f64(dotproduct);
//#else
        return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
//#endif        
    }
    inline static double3 max(const double3& v1, const double3& v2) { return double3(std::max(v1.x, v2.x), std::max(v1.y, v2.y), std::max(v1.z, v2.z)); }
    inline static double3 min(const double3& v1, const double3& v2) { return double3(std::min(v1.x, v2.x), std::min(v1.y, v2.y), std::min(v1.z, v2.z)); }
    inline static double3 cross(const double3& v1, const double3& v2) { return double3(v1.y * v2.z - v2.y * v1.z, v1.z * v2.x - v2.z * v1.x, v1.x * v2.y - v2.x * v1.y); }
    inline static double3 normalize(const double3& v) { double f = 1.0 / sqrt((v.x * v.x) + (v.y * v.y) + (v.z * v.z)); return double3(f * v.x, f * v.y, f * v.z); }
    inline static  double3 flip(double3 v, bool3 flip, double3 boundary) {
        return double3(flip.x ? boundary.x - v.x : v.x,
            flip.y ? boundary.y - v.y : v.y,
            flip.z ? boundary.z - v.z : v.z);
    }
    

    double3 operator-() const { return double3(-this->x, -this->y, -this->z); }
    double3& operator+=(const double3& b) { this->x += b.x, this->y += b.y, this->z += b.z; return *this; }
    double3& operator-=(const double3& b) { this->x -= b.x, this->y -= b.y, this->z -= b.z; return *this; }

    friend std::ostream& operator<<(std::ostream& out, const double3& vec);

    struct KeyHash {
        size_t operator()(const double3& k) const {
            size_t h1 = std::hash<double>()(k.x);
            size_t h2 = std::hash<double>()(k.y);
            size_t h3 = std::hash<double>()(k.z);
            return (h1 ^ (h2 << 1)) ^ h3;
        }
    };

    struct KeyEqual {
        bool operator()(const double3& lhs, const double3& rhs) const {
            return (std::fabs(lhs.x - rhs.x) < 1e-3) && std::fabs(lhs.y - rhs.y) < 1e-3 && std::fabs(lhs.z == rhs.z) < 1e-3;
        }
    };
};


export inline double3 operator+(const double3& a, const double3& b)
{
//#ifdef __AVX__
//    return double3(_mm256_add_pd(a.vec, b.vec));
//#else
    return double3(a.x + b.x, a.y + b.y, a.z + b.z);
//#endif
}

export inline double3 operator-(const double3& a, const double3& b)
{
//#ifdef __AVX__
//    return double3(_mm256_sub_pd(a.vec, b.vec));
//#else
    return double3(a.x - b.x, a.y - b.y, a.z - b.z);
//#endif
}

export inline double3 operator*(const double3& a, const double3& b)
{
    return double3(a.x * b.x, a.y * b.y, a.z * b.z);
}

export inline double3 operator/(const double3& a, const double3& b)
{
    return double3(a.x / b.x, a.y / b.y, a.z / b.z);
}

export inline double3 operator*(const double3& a, double b)
{
    return double3(a.x * b, a.y * b, a.z * b);
}

export inline double3 operator*(const double& a, const double3& b)
{
    return double3(a * b.x, a * b.y, a * b.z);
}

export inline double3 operator/(const double3& a, double b)
{
    return double3(a.x / b, a.y / b, a.z / b);
}

export inline double3 sqrt(const double3& a)
{
    return double3(std::sqrt(a.x), std::sqrt(a.y), std::sqrt(a.z));
}
