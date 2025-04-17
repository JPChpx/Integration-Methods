 # 🧮 Numerical Integration Library

A lightweight and powerful C++ library for numerical integration methods. Developed by JP Champeaux (2019), this library includes a variety of 1D, 2D, and nD integration algorithms — from classical rules like Trapezoidal and Simpson, to advanced techniques like Romberg and Monte Carlo.

---

## 📦 Features

- ✅ Rectangular (left, middle, right)
- ✅ Trapezoidal Rule
- ✅ Simpson’s 1/3 Rule
- ✅ Simpson’s 3/8 Rule
- ✅ Bode’s Rule
- ✅ Weddle’s Rule
- ✅ High-order Rules (7th to 10th order)
- ✅ Romberg Integration
- ✅ 2D Simpson Integration (data + continuous functions)
- ✅ nD Monte Carlo Integration

---

## 📁 File Structure

- `integration.h` – All implementations are included in this header.
- *(Optional)* Add a `main.cpp` file to test and benchmark methods.

---

## ⚙️ Usage

###  Example – Using Simpson’s Rule:
```cpp
auto f = [](double x) { return sin(x); };
double result = Std_integrate(simpson(), f, 0.0, M_PI, 1000);
std::cout << "Integral of sin(x) from 0 to PI ≈ " << result << std::endl;
```

###  Example – Monte Carlo nD Integration:
```cpp
auto g = [](std::vector<double> x) { return x[0] * x[1]; };
std::vector<double> a = {0, 0}, b = {1, 1};
double result = MC_integrate(g, a, b, 1e-4, 100000);
```
---
# 🛠 Requirements

    C++11 or later
    Standard C++ libraries (no external dependencies)
---
#📜 License

MIT License
© JP Champeaux 2019 – Modified and structured for open-source sharing.
🤝 Contributions

Pull requests, bug reports, or feature suggestions are welcome!
Let's improve this library together 🚀
