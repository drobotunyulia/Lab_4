#include "longint.h"

long_int::long_int()
{
}

long_int::long_int(int value)
{
    this->sign = value < 0;
    std::string string_value = std::to_string(value);
    for(size_t i = 0 + static_cast<size_t>(sign); i < string_value.length(); i++)
    {
        this->digits.push_back(string_value[i] - '0');
    }
}

long_int::long_int(const std::string &value, size_t radix)
{
    if((radix < 2) || (radix > 36))
    {
        throw "Incorrect radix of number systems.";
    }
    this->sign = value[0] == '-';
    if(radix == 10)
    {
        for(size_t i = 0 + static_cast<size_t>(sign); i < value.length(); i++)
        {
            if(std::isdigit(value[i]))
            {
                this->digits.push_back(value[i] - '0');
            }
            else
            {
                throw "Invalid character when writing a number.";
            }
        }
    }
    else
    {
        long_int raw_value;
        long_int long_radix = radix;
        long_int result;
        long_int digit;
        long_int degree_factor;
        for(size_t i = 0 + static_cast<size_t>(sign); i < value.length(); i++)
        {
            int digit;
            if(std::isdigit(value[i]))
            {
                digit = value[i] - 48;
            }
            else if(std::isupper(value[i]) && radix > 10)
            {
                digit = value[i] - 55;
            }
            else if(std::islower(value[i]) && radix > 10)
            {
                digit = value[i] - 87;
            }
            else
            {
                throw "Invalid character when writing a number.";
            }
            raw_value.digits.push_back(digit);
        }
        while(raw_value.digits[0] == 0 && raw_value.size() > 1)
        {
            raw_value.digits.erase(raw_value.digits.begin());
        }
        std::reverse(raw_value.digits.begin(), raw_value.digits.end());
        if(raw_value != 0)
        {
            for(size_t i = 0; i < raw_value.size(); i++)
            {
                degree_factor = long_radix.fast_pow(i);
                digit = raw_value.digits[i];
                digit = digit.mul(degree_factor);
                result = result.add(digit);
            }
            this->digits = result.digits;
        }
        else
        {
            *this = 0;
        }
    }
}

long_int::long_int(const long_int &long_int_obj)
{
    this->digits = long_int_obj.digits;
    this->sign = long_int_obj.sign;
}

long_int::~long_int()
{
    this->digits.clear();
}

long_int long_int::add(const long_int &value) const
{
    if(*this == 0 && value == 0)
    {
        return 0;
    }
    long_int result;
    if (!(value.sign ^ this->sign))
    {
        std::vector<digit_t> digits_2 = value.digits;
        size_t size_digits_1 = this->size();
        size_t size_digits_2 = digits_2.size();
        size_t size_result = 1 + std::max(size_digits_1, size_digits_2);
        result.digits.resize(size_result);
        for (size_t i = 0; i < size_result - 1; i++)
        {
            int j = size_result - 1 - i;
            result.digits[j] += ((i < size_digits_2) ? (digits_2[size_digits_2 - 1 - i]) : 0) +
                    ((i < size_digits_1) ? (this->digits[size_digits_1 - 1 - i]) : 0);
            result.digits[j - 1] = result.digits[j] / 10;
            result.digits[j] = result.digits[j] % 10;
        }
        result.sign = this->sign;
        while(result.digits[0] == 0 && result.size() > 1)
        {
            result.digits.erase(result.digits.begin());
        }
    }
    else
    {
        result = this->sign ? (value.sub(-long_int(*this))) : (this->sub(-long_int(value)));
    }
    return result;
}

long_int long_int::sub(const long_int &value) const
{
    long_int result;
    if (*this == value)
    {
        return 0;
    }
    if (!this->sign && !value.sign)
    {
        std::vector<digit_t> digits_2 = value.digits;
        size_t size_digits_1 = this->size();
        size_t size_digits_2 = digits_2.size();
        size_t size_result = std::max(size_digits_1, size_digits_2);
        result.digits.clear();
        result.digits.resize(size_result);
        bool sign_result = value > *this;
        std::vector<digit_t> inner_value_1(size_result), inner_value_2(size_result);
        inner_value_1[0] = inner_value_2[0] = 0;
        int sign = (2 * sign_result - 1);
        for (size_t i = 0; i < size_result - 1; i++)
        {
            inner_value_1[i] += (i < size_digits_1) ? (this->digits[size_digits_1 - 1 - i]) : 0;
            inner_value_2[i] += (i < size_digits_2) ? (digits_2[size_digits_2 - 1 - i]) : 0;
            inner_value_2[i + 1] = -sign_result;
            inner_value_1[i + 1] = sign_result - 1;
            result.digits[size_result - 1 - i] += 10 + sign * (inner_value_2[i] - inner_value_1[i]);
            result.digits[size_result - 1 - i - 1]  = result.digits[size_result - 1 - i] / 10;
            result.digits[size_result - 1 - i] = result.digits[size_result - 1 - i] % 10;
        }
        inner_value_1[size_result - 1] += (size_result - 1 < size_digits_1) * (this->digits[0]);
        inner_value_2[size_result - 1] += (size_result - 1 < size_digits_2) * (digits_2[0]);
        result.digits[0] += sign * (inner_value_2[size_result - 1] - inner_value_1[size_result - 1]);
        result.sign = sign_result;
        while(result.digits[0] == 0)
        {
            result.digits.erase(result.digits.begin());
        }
    }
    else
    {
        result = this->sign && value.sign ? (-long_int(value).sub(-long_int(*this))) : (this->add(-long_int(value)));
        while(result.digits[0] == 0 && result.size() > 1)
        {
            result.digits.erase(result.digits.begin());
        }
    }
    return result;
}

long_int long_int::mul(const long_int &value) const
{
    if(*this == 0 || value == 0)
    {
        return 0;
    }
    std::vector<digit_t> digits_2 = value.digits;
    size_t size_digits_1 = this->size();
    size_t size_digits_2 = digits_2.size();
    size_t size_result = size_digits_1 + size_digits_2 + 1;
    bool sign_result = this->sign ^ value.sign;
    std::vector<digit_t> inner_value_1(size_result), inner_value_2(size_result);
    long_int result;
    result.digits.resize(size_result);
    for (size_t i = 0; i < size_result; i++)
    {
        inner_value_1[i] = (i < size_digits_1) ? (this->digits[size_digits_1 - 1 - i]) : 0;
        inner_value_2[i] = (i < size_digits_2) ? (digits_2[size_digits_2 - 1 - i]) : 0;
    }
    for (size_t i = 0; i < size_digits_1; i++)
    {
        for (size_t j = 0; j < size_digits_2; j++)
        {
            result.digits[size_result - 1 - (i + j)] += inner_value_1[i] * inner_value_2[j];
            result.digits[size_result - 1 - (i + j + 1)] += result.digits[size_result - 1 - (i + j)] / 10;
            result.digits[size_result - 1 - (i + j)] %= 10;
        }
    }
    result.sign = sign_result;
    while(result.digits[0] == 0 && result.size() > 1)
    {
        result.digits.erase(result.digits.begin());
    }
    return result;
}

long_int long_int::mul_karatsuba(const long_int &value) const
{
    bool sign_result = this->sign ^ value.sign;
    std::vector<digit_t> digits_1 = this->digits;
    std::vector<digit_t> digits_2 = value.digits;
    long_int value_1 = *this;
    long_int value_2 = value;
    size_t size_digits_1 = digits_1.size();
    size_t size_digits_2 = digits_2.size();
    size_t size_max = std::max(size_digits_1, size_digits_2);
    if (size_digits_1 == 1 || size_digits_2 == 1)
    {
        return long_int(value_1.mul(value_2));
    }
    size_max += size_max % 2;
    size_t size_half = size_max / 2;
    long_int value_1_right;
    size_digits_1 > size_half ? value_1_right.digits.assign(digits_1.end() - size_half, digits_1.end()) :
                        value_1_right.digits.assign(digits_1.begin(), digits_1.end());
    long_int value_2_right;
    size_digits_2 > size_half ? value_2_right.digits.assign(digits_2.end() - size_half, digits_2.end()) :
                        value_2_right.digits.assign(digits_2.begin(), digits_2.end());
    long_int value_1_left(value_1.shr(size_half));
    long_int value_2_left(value_2.shr(size_half));
    long_int P1 = value_1_left.mul_karatsuba(value_2_left);
    long_int P2 = value_1_right.mul_karatsuba(value_2_right);
    long_int P3 = value_1_left.add(value_1_right).mul_karatsuba(value_2_left.add(value_2_right));
    long_int result = P1.shl(size_max).add((P3.sub(P2).sub(P1).shl(size_half)).add(P2));
    result.sign = sign_result;
    return result;
}

long_int long_int::mul_fft(const long_int &value) const
{
    bool sign_result = this->sign ^ value.sign;
    long_int result;
    std::vector<std::complex<double>> complex_value_1(this->digits.begin(), this->digits.end()),
            complex_value_2(value.digits.begin(), value.digits.end());
    size_t size_result = 1;
    std::reverse(complex_value_1.begin(), complex_value_1.end());
    std::reverse(complex_value_2.begin(), complex_value_2.end());
    while (size_result < std::max(complex_value_1.size(), complex_value_2.size()))
    {
        size_result <<= 1;
    }
    size_result <<= 1;
    complex_value_1.resize(size_result);
    complex_value_2.resize(size_result);

    complex_value_1 = this->get_fft(complex_value_1, false);
    complex_value_2 = this->get_fft(complex_value_2, false);
    for (size_t i = 0; i < size_result; ++i)
    {
        complex_value_1[i] *= complex_value_2[i];
    }
    complex_value_1 = this->get_fft (complex_value_1, true);
    int carry = 0;
    digit_t result_digit;
    for (size_t i = 0; i < size_result; ++i)
    {
        result_digit = static_cast<digit_t>(complex_value_1[i].real() + 0.5) + carry;
        carry = result_digit / 10;
        result_digit %= 10;
        result.digits.push_back(result_digit);
    }
    result.digits.push_back(carry);
    std::reverse(result.digits.begin(), result.digits.end());
    while(result.digits[0] == 0 && result.size() > 1)
    {
        result.digits.erase(result.digits.begin());
    }
    result.sign = sign_result;
    return result;
}

long_int long_int::div_newton_raphson(const long_int &value) const
{
    bool sign_result = this->sign ^ value.sign;
    if(value == 0)
    {
        throw "Division by zero.";
    }
    if(value > *this)
    {
        return 0;
    }
    if(value == *this)
    {
        return 1;
    }
    if(value == 1)
    {
        return *this;
    }
    long_int divident = *this;
    long_int divider = value;
    size_t shift = (divident.size() + divider.size());
    long_int pow_10(2);
    pow_10 = pow_10.shl(shift);
    long_int x_current = divident.sub(divider);
    long_int x_last;
    while(x_last != x_current)
    {
        x_last = x_current;
        x_current = (x_current.mul(pow_10.sub(x_current.mul(divider)))).shr(shift);
    }
    long_int result = divident.mul(x_current);
    result = result.shr(shift);
    long_int remainder = divident.sub(result.mul(divider));
    if(remainder >= divider)
    {
        result = result.add(1);
    }
    result.sign = sign_result;
    return result;
}

long_int long_int::fast_pow(size_t degree) const
{
    if (degree == 0)
    {
        return 1;
    }
    long_int inner_value;
    if (degree & 1)//провекра на нечетность
    {
        inner_value = this->fast_pow(degree - 1);
        return inner_value.mul(*this);
    }
    else
    {
        inner_value = this->fast_pow(degree / 2);
        return inner_value.mul(inner_value);
    }
}

long_int long_int::fast_pow_mod(size_t degree, const long_int mod) const
{
    if (degree == 0)
    {
        return 1;
    }
    long_int inner_value;
    if (degree & 1)
    {
        inner_value = this->fast_pow_mod(degree - 1, mod);
        return inner_value.mul(*this).rem(mod);
    }
    else
    {
        inner_value = this->fast_pow_mod(degree / 2, mod);
        return inner_value.mul(inner_value).rem(mod);
    }
}

long_int long_int::rem(const long_int &value) const
{
    if(value > *this)
    {
        return *this;
    }
    if(value == *this)
    {
        return 0;
    }
    long_int divident(*this);
    long_int divider(value);
    long_int quotient = divident.div_newton_raphson(divider);
    return divident.sub(divider.mul(quotient));
}

long_int long_int::shl(size_t sh_num) const
{
    long_int result(*this);
    for(size_t i = 0; i < sh_num; i++)
    {
        result.digits.push_back(0);
    }
    return result;
}

long_int long_int::shr(size_t sh_num) const
{
    if(sh_num >= this->size())
    {
        return long_int(0);
    }
    long_int result(*this);
    for(size_t i = 0; i < sh_num; i++)
    {
        result.digits.erase(result.digits.end() - 1);
    }
    return result;
}

long_int& long_int::operator = (const long_int &long_int_obj)
{
    this->digits = long_int_obj.digits;
    this->sign = long_int_obj.sign;
    return *this;
}

long_int& long_int::operator = (int value)
{
    *this = long_int(value);
    return *this;
}

long_int long_int::operator-() const &
{
    long_int result(*this);
    if(this->sign)
    {
        result.sign = false;
    }
    else
    {
        result.sign = true;
    }
    return result;
}

long_int long_int::operator+() const &
{
    return long_int(*this);
}

bool long_int::operator == (const long_int &long_int_obj) const
{
    return (this->digits == long_int_obj.digits) && (sign == long_int_obj.sign);
}

bool long_int::operator != (const long_int &long_int_obj) const
{
    return !(*this == long_int_obj);
}

bool long_int::operator < (const long_int &long_int_obj) const
{
    std::vector<digit_t> digits_2 = long_int_obj.digits;
    size_t size_digits_1 = this->size();
    size_t size_digits_2 = digits_2.size();

    if (this->sign == long_int_obj.sign)
    {
        if (size_digits_1 != size_digits_2)
            return (size_digits_1 < size_digits_2) ^ this->sign;
        size_t i = 0;
        while (i < size_digits_1 && this->digits[i] == digits_2[i])
        {
            i++;
        }
        return (i < size_digits_1) && ((this->digits[i] < digits_2[i]) ^ this->sign);
    }
    return this->sign;
}

bool long_int::operator > (const long_int &long_int_obj) const
{
    return !(*this < long_int_obj || *this == long_int_obj);
}

bool long_int::operator <= (const long_int &long_int_obj) const
{
    return *this < long_int_obj || *this == long_int_obj;
}

bool long_int::operator >= (const long_int &long_int_obj) const
{
    return *this > long_int_obj || *this == long_int_obj;
}

bool long_int::is_even() const
{
    if(this->digits[this->size() - 1] & 1)
    {
        return false;
    }
    return true;
}

size_t long_int::num_bits() const
{
    long_int inner_value(*this);
    size_t result = 0;
    if(*this == 0)
    {
        return 0;
    }
    while(inner_value != 0)
    {
        result++;
        inner_value = inner_value.div_newton_raphson(2);
    }
    return result;
}


size_t long_int::size() const
{
    return this->digits.size();
}

std::string long_int::to_string() const
{
    std::string result;
    if(sign)
    {
        result.append("-");
    }
    for(size_t i = 0; i < this->size(); i++)
    {
        result.append(std::to_string(this->digits[i]));
    }
    return result;
}

std::vector<std::complex<double>> long_int::get_fft(const std::vector<std::complex<double>> &value,
                                                    bool invert) const
{
    std::vector<std::complex<double>> result = value;
    int size_result = (int)value.size();
    if (size_result == 1)
    {
        return std::vector<std::complex<double>>(1, value[0]);
    }
    std::vector<std::complex<double>> value_not_even(size_result / 2);
    std::vector<std::complex<double>> value_even(size_result / 2);
    for (int i = 0, j = 0; i < size_result; i += 2, ++j)
    {
        value_not_even[j] = value[i];
        value_even[j] = value[i+1];
    }
    value_not_even = get_fft(value_not_even, invert);
    value_even = get_fft(value_even, invert);
    double angle = 2 * PI / size_result * (invert ? -1 : 1);
    std::complex<double> root_1(1);
    std::complex<double> root_n(std::cos(angle), std::sin(angle));
    for (int i = 0; i < size_result/2; ++i)
    {
        result[i] = value_not_even[i] + root_1 * value_even[i];
        result[i + size_result / 2] = value_not_even[i] - root_1 * value_even[i];
        if (invert)
        {
            result[i] /= 2;
            result[i + size_result / 2] /= 2;
        }
        root_1 *= root_n;
    }
    return result;
}

long_int gcd(const long_int &value_1, const long_int &value_2)
{
    long_int inner_value_1 = std::max(value_1, value_2);
    long_int inner_value_2 = std::min(value_1, value_2);
    while (inner_value_2 > 0)
    {
        inner_value_1 = inner_value_1.rem(inner_value_2);
        std::swap(inner_value_1, inner_value_2);
    }
    return inner_value_1;

}

long_int gcd_ext(const long_int &value_1, const long_int &value_2,
                        long_int &factor_1, long_int &factor_2)
{
    long_int inner_value_1 = std::max(value_1, value_2);
    long_int inner_value_2 = std::min(value_1, value_2);
    long_int q, r, x_1, x_2, y_1, y_2, result;
    if (inner_value_2 == 0)
    {
        result = inner_value_1;
        factor_1 = 1;
        factor_2 = 0;
        return result;
    }
    x_2 = 1, x_1 = 0, y_2 = 0, y_1 = 1;
    while (inner_value_2 > 0)
    {
        q = inner_value_1.mul(inner_value_2);
        r = inner_value_1.sub(q.mul(inner_value_2));
        factor_1 = x_2.sub(q.mul(x_1));
        factor_2 = y_2.sub(q.mul(y_1));
        inner_value_1 = inner_value_2, inner_value_2 = r;
        x_2 = x_1;
        x_1 = factor_1;
        y_2 = y_1;
        y_1 = factor_2;
    }
    result = inner_value_1;
    factor_1 = x_2;
    factor_2 = y_2;
    return result;
}

long_int gcd_binary(const long_int &value_1, const long_int &value_2)
{
    long_int inner_value_1 = std::max(value_1, value_2);
    long_int inner_value_2 = std::min(value_1, value_2);
    long_int k = 1;
    while ((inner_value_1 != 0) && (inner_value_2 != 0))
    {
        while ((inner_value_1.is_even()) && (inner_value_2.is_even()))
        {
            inner_value_1 = inner_value_1.div_newton_raphson(2);
            inner_value_2 = inner_value_2.div_newton_raphson(2);
            k = k.mul(2);
        }
        while (inner_value_1.is_even())
        {
            inner_value_1 = inner_value_1.div_newton_raphson(2);
        }
        while (inner_value_2.is_even())
        {
            inner_value_2 = inner_value_2.div_newton_raphson(2);
        }
        if (inner_value_1 >= inner_value_2)
        {
            inner_value_1 = inner_value_1.sub(inner_value_2);
        }
        else
        {
            inner_value_2 = inner_value_2.sub(inner_value_1);
        }
    }
    return inner_value_2.mul(k);
}

int get_legendre_symbol(const long_int &value_a, const long_int &value_p)
{
    long_int a_mod_p = value_a.rem(value_p);
    if (a_mod_p == 0)
    {
        return 0;
    }
    if (a_mod_p == 1)
    {
        return 1;
    }
    for (long_int i = 1; i <= value_p; i = i.add(1))
    {
        if (i.mul(i).rem(value_p) == a_mod_p)
        {
            return 1;
        }
    }
    return -1;
}

int get_jacobi_symbol(const long_int &value_a, const long_int &value_p)
{
    if(value_p.is_even() || value_p < 1)
    {
        throw "Invalid value of function parameters.";
    }
    if(value_p == 1)
    {
        return 1;
    }
    long_int inner_value_a = value_a;
    long_int inner_value_p = value_p;
    if (gcd(inner_value_a, inner_value_p) != 1)
    {
        return 0;
    }
    int result = 1;

    if (inner_value_a < 0)
    {
        inner_value_a = -inner_value_a;
        if ((inner_value_p.rem(4)) == 3)
        {
                result = -result;
        }
    }
    while (inner_value_a != 0)
    {
        long_int temp_value_1 = 0;
        while (inner_value_a.is_even())
        {
            temp_value_1 = temp_value_1.add(1);
            inner_value_a = inner_value_a.div_newton_raphson(2);
        }
        if (!temp_value_1.is_even())
        {
            if (((inner_value_p.rem(8)) == 3) || ((inner_value_p.rem(8)) == 5))
            {
                result = -result;
            }
        }
        if ((inner_value_a.rem(4) == inner_value_p.rem(4)) && (inner_value_p.rem(4) == 3))
        {
            result = -result;
        }
        long_int temp_value_2 = inner_value_a;
        inner_value_a = inner_value_p.rem(temp_value_2);
        inner_value_p = temp_value_2;
    }
    return result;
}

bool check_primary_fermat(const long_int &value_value_a, const long_int &value_p)
{

}

bool check_primary_solovay_strassen(const long_int &value_a, long_int &value_p)
{

}

bool check_primary_miller_rabin(const long_int &value_a, long_int &value_p)
{

}

std::ostream& operator << (std::ostream &out, const long_int &long_int_obj)
{
    out << long_int_obj.to_string();
    return out;
}

rand_generator::rand_generator()
{
    srand(time(0));
}

long_int *rand_generator::generate(size_t num_bits) const
{
    long_int raw_value;
    long_int *result = new long_int(0);
    long_int hex_radix(16);
    long_int digit;
    long_int degree_factor;
    size_t remainder = num_bits % 8;

    unsigned char bit_mask;
    for(size_t i = 0; i < num_bits / 4; i++)
    {
        raw_value.digits.push_back(this->get_random_digit());
    }
    if(remainder != 0)
    {
        bit_mask = static_cast<unsigned char>(std::pow(2, remainder) - 1);
        raw_value.digits.push_back(this->get_random_digit() & bit_mask);
    }
    std::reverse(raw_value.digits.begin(), raw_value.digits.end());
    if(raw_value != 0)
    {
        for(size_t i = 0; i < raw_value.size(); i++)
        {
            degree_factor = hex_radix.fast_pow(i);
            digit = raw_value.digits[i];
            digit = digit.mul(degree_factor);
            *result = (*result).add(digit);
        }
    }
    return result;
}

size_t rand_generator::get_random_digit() const
{
    static const double fraction = 1.0 / (static_cast<double>(RAND_MAX) + 1.0);
    return static_cast<size_t>(rand() * fraction * (16));
}
