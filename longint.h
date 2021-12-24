#ifndef LONGINT_H
#define LONGINT_H

#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cctype>
#include <iostream>
#include <vector>
#include <complex>
#include <string>
#include <algorithm>

using digit_t = short int;

const double PI = std::acos(-1);

class rand_generator;

class long_int
{
private:
    std::vector<digit_t> digits;
    bool sign = false;

    std::vector<std::complex<double>> get_fft(const std::vector<std::complex<double> > &value,
                                              bool invert) const;
public:
    long_int();
    long_int(int value);
    long_int(const std::string &value, size_t radix = 10);//конструктор от строки с числом в заданной системе счисления
    long_int(const long_int &long_int_obj);//конструктор копирования
    ~long_int();

    long_int add(const long_int &value) const;//сложение чисел
    long_int sub(const long_int &value) const;//вычитание чисел
    long_int mul(const long_int &value) const;//умножение в столбик
    long_int mul_karatsuba(const long_int &value) const;//алгоритм Карацубы
    long_int mul_fft(const long_int &value) const;//умножение чисел (на основе быстрого преобразования Фурье)
    long_int div_newton_raphson(const long_int &value) const;//метод ньютона
    long_int fast_pow(size_t degree) const;//быстрое возведение в степень
    long_int fast_pow_mod(size_t degree, const long_int mod) const;//быстрое возведение в степень по модулю
    long_int rem(const long_int &value) const;//остаток

    long_int shl(size_t sh_num) const;//сдвиг влево
    long_int shr(size_t sh_num) const;//сдвиг вправо

    long_int& operator = (const long_int &long_int_obj);
    long_int& operator = (int value);

    long_int operator+() const &;
    long_int operator-() const &;

    bool operator == (const long_int &long_int_obj) const;
    bool operator != (const long_int &long_int_obj) const;
    bool operator < (const long_int &long_int_obj) const;
    bool operator > (const long_int &long_int_obj) const;
    bool operator <= (const long_int &long_int_obj) const;
    bool operator >= (const long_int &long_int_obj) const;

    size_t num_bits() const;//количество бит
    bool is_even() const;//проверка на четность
    size_t size() const;//размер числа
    std::string to_string() const;//преобразование в строку

    friend long_int gcd(const long_int &value_1, const long_int &value_2);
    friend long_int gcd_ext(const long_int &value_1, const long_int &value_2,
                            long_int &factor_1, long_int &factor_2);
    friend long_int gcd_binary(const long_int &value_1, const long_int &value_2);

    friend int get_legendre_symbol(const long_int &value_a, const long_int &value_p);
    friend int get_jacobi_symbol(const long_int &value_a, const long_int &value_p);
    friend bool check_primary_fermat(const long_int &value_a, const long_int &value_p);
    friend bool check_primary_solovay_strassen(const long_int &value_a, long_int &value_p);
    friend bool check_primary_miller_rabin(const long_int &value_a, long_int &value_p);

    friend std::ostream& operator << (std::ostream &out, const long_int &long_int_obj);
    friend class rand_generator;
};

long_int gcd(const long_int &value_1, const long_int &value_2);
long_int gcd_ext(const long_int &value_1, const long_int &value_2,
                 long_int &factor_1, long_int &factor_2);
long_int gcd_binary(const long_int &value_1, const long_int &value_2);

int get_legendre_symbol(const long_int &value_a, const long_int &value_p);//вычисление символа Лежандра
int get_jacobi_symbol(const long_int &value_a, const long_int &value_p);//вычисление символа Якоби
bool check_primary_fermat(const long_int &value_a, const long_int &value_p);//определение простоты числа при помощи теста Ферма
bool check_primary_solovay_strassen(const long_int &value_a, long_int &value_p);//определение простоты числа при помощи теста Соловея-Штрассена
bool check_primary_miller_rabin(const long_int &value_a, long_int &value_p);//определение простоты числа при помощи теста Миллера-Рабина

std::ostream& operator << (std::ostream &out, const long_int &long_int_obj);


class rand_generator
{
private:
    size_t get_random_digit() const;
public:
    rand_generator();
    long_int *generate(size_t num_bits) const;
};

using long_int_t = long_int;
using rand_generator_t = rand_generator;

#endif // LONGINT_H
