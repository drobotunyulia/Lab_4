#include <iostream>
#include <cmath>

#include "longint.h"

using namespace std;

int main()
{
    long_int_t value_1("9000000000000");
    cout << "Val 1: " << value_1 << endl;
    long_int_t value_2("4000fff0", 16);
    cout << "Val 2: " << value_2 << endl;
    long_int_t value_3("3456436234623");
    cout << "Val 3: " << value_3 << endl;
    long_int_t value_4("3426534667300");
    cout << "Val 4: " << value_4 << endl;
    long_int_t value_5("3456436234623");
    cout << "Val 5: " << value_5 << endl;
    long_int_t value_6("3426534667300");
    cout << "Val 6: " << value_6 << endl;
    cout << "Val 1 * Val 2 = " << value_1.mul_fft(value_2) << endl;
    cout << " Val 3 * Val 4 = " << value_3.mul(value_4) << endl;
    cout << "Val 5 * Val 6 = "<< value_5.mul_karatsuba(value_6) << endl;
    cout << "Val 1 - Val 2 = " << value_1.sub(value_2) << endl;
    cout << "Val 1 / 666666666 = " << value_1.div_newton_raphson(666666666) << endl;
    long_int_t value_7("595789");
    cout << "Val 7: " << value_7 << endl;
    long_int_t value_8("595789");
    cout << "Val 8: " << value_8 << endl;
    cout << "Reminder = " << value_7.rem(value_8) << endl;
    rand_generator_t *rand_generator = new rand_generator_t;
    cout << *rand_generator->generate(256) << endl;
    cout << *rand_generator->generate(256) << endl;
    cout << *rand_generator->generate(256) << endl;
    cout << (*rand_generator->generate(96)).num_bits() << endl;
    cout <<  "GCD Val 1 and Val 2 = " <<gcd(value_1, value_2) << endl;
    cout << gcd_binary(value_1, value_2) << endl;
    cout << get_legendre_symbol(456, 79) << endl;
    return 0;
}
