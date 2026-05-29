import unittest
from factorize_by_digits import sieve, m_digit_primes, factorize_by_digits


class TestSieve(unittest.TestCase):
    def test_sieve_small(self):
        self.assertEqual(sieve(10), [2, 3, 5, 7])

    def test_sieve_30(self):
        self.assertEqual(sieve(30), [2, 3, 5, 7, 11, 13, 17, 19, 23, 29])

    def test_sieve_2(self):
        self.assertEqual(sieve(2), [2])


class TestMDigitPrimes(unittest.TestCase):
    def test_1_digit(self):
        self.assertEqual(m_digit_primes(1), [2, 3, 5, 7])

    def test_2_digit(self):
        res = m_digit_primes(2)
        self.assertEqual(len(res), 21)
        self.assertEqual(res[0], 11)
        self.assertEqual(res[-1], 97)


class TestFactorizeByDigits(unittest.TestCase):
    def test_11x13(self):
        self.assertEqual(factorize_by_digits(143, 2), [(11, 13)])

    def test_17x19(self):
        self.assertEqual(factorize_by_digits(323, 2), [(17, 19)])

    def test_29x31(self):
        self.assertEqual(factorize_by_digits(899, 2), [(29, 31)])

    def test_31x37(self):
        self.assertEqual(factorize_by_digits(1147, 2), [(31, 37)])

    def test_41x43(self):
        self.assertEqual(factorize_by_digits(1763, 2), [(41, 43)])

    def test_11x11(self):
        self.assertEqual(factorize_by_digits(121, 2), [(11, 11)])

    def test_no_factors_prime(self):
        self.assertEqual(factorize_by_digits(211, 2), [])


if __name__ == "__main__":
    unittest.main()
