#include <iostream>
#include <iomanip>
#include <cmath>
#include <functional>
#include <string>
#include <array>
using namespace std;

using Func = function<long double(long double)>;

class NewCotes {
public:
    NewCotes(Func f) : f(f) {}

    long double integrale(long double a, long double b, int n) const {
        long double h = (b - a) / n;
        long double s = 0;
        for (int i = 0; i < n; ++i) {
            long double x = a + i * h;
            s += h / 6 * (f(x) + 4 * f(x + h / 2) + f(x + h));
        }
        return s;
    }

    int order() const { return 4; }

private:
    Func f;
};

class Gauss {
public:
    Gauss(Func f) : f(f) {}

    long double integrale(long double a, long double b, int n) const {
        long double h = (b - a) / n;
        long double d = h / (2 * sqrtl(3.0L));   
        long double s = 0;
        for (int i = 0; i < n; ++i) {
            long double c = a + i * h + h / 2;  
            s += h / 2 * (f(c - d) + f(c + d));
        }
        return s;
    }

    int order() const { return 4; }

private:
    Func f;
};

Func punkt2_funkciya() {
    return [](long double x) { return expl(x); };      // phi(x) = e^x
}

long double punkt3_tochnoe(long double a, long double b) {
    return expl(b) - expl(a);                          // e^b - e^a
}



template <class Scheme>
array<long double, 3> integralyNaSetkah(const Scheme& q, long double a, long double b, int n0) {
    return { q.integrale(a, b, n0),
             q.integrale(a, b, 2 * n0),
             q.integrale(a, b, 4 * n0) };
}



void punkt4_setki(long double a, long double b, int n0) {
    cout << "Setki: h = " << (b - a) / n0
         << ", h/2 = " << (b - a) / (2 * n0)
         << ", h/4 = " << (b - a) / (4 * n0) << "\n\n";
}



template <class Scheme>
void punkt4_poryadok(const string& name, const Scheme& q, long double a, long double b, int n0) {
    auto I = integralyNaSetkah(q, a, b, n0);
    long double k = log2l((I[0] - I[1]) / (I[1] - I[2]));
    cout << name << ": I(h) = " << I[0] << ", I(h/2) = " << I[1]
         << ", I(h/4) = " << I[2] << "\n";
    cout << name << ": k = " << k << "\n\n";
}



template <class Scheme>
void punkt5_tablica(const string& name, const Scheme& q, long double a, long double b,
                    int n0, long double Istar) {
    auto I = integralyNaSetkah(q, a, b, n0);
    int k = q.order();
    long double c = powl(2.0L, k) - 1;                 
    long double h0 = (b - a) / n0;

    cout << name << " (k = " << k << ")\n";
    cout << "h | I*-Ih | (I*-Ih)/(I*-Ih/2) | (Ih/2-Ih)/(2^k-1) | IR | I*-IR\n";
    for (int j = 0; j < 2; ++j) {                     
        long double h     = h0 / powl(2.0L, j);
        long double err   = Istar - I[j];              // I* - I^h
        long double ratio = err / (Istar - I[j + 1]);  // (I*-I^h)/(I*-I^{h/2})
        long double runge = (I[j + 1] - I[j]) / c;     // (I^{h/2}-I^h)/(2^k-1)
        long double IR    = I[j + 1] + runge;          
        cout << h << "  " << err << "  " << ratio << "  "
             << runge << "  " << IR << "  " << Istar - IR << "\n";
    }
    cout << "\n";
}


int main() {
    const long double a = 0.03L, b = 0.30L;            
    const int n0 = 4;                                  
    cout << scientific << setprecision(12);

    Func phi = punkt2_funkciya();                      
    long double Istar = punkt3_tochnoe(a, b);          
    cout << "I* = " << Istar << "\n\n";

    NewCotes simpson(phi);
    Gauss gauss(phi);

    punkt4_setki(a, b, n0);                          
    punkt4_poryadok("Simpson", simpson, a, b, n0);
    punkt4_poryadok("Gauss2", gauss, a, b, n0);

    punkt5_tablica("Simpson", simpson, a, b, n0, Istar);   
    punkt5_tablica("Gauss2", gauss, a, b, n0, Istar);
    return 0;
}