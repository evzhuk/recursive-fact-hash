import unittest
from factorize_by_digits import sieve, m_digit_primes, factorize_by_digits
from prime_utils import (
    passes_small_prime_filter,
    is_probable_prime_mr,
    is_prime_fast,
    m_digit_primes_fast,
    PrimeStats,
)


class TestSieve(unittest.TestCase):
    """Legacy sieve (не изменился)."""

    def test_sieve_small(self):
        self.assertEqual(sieve(10), [2, 3, 5, 7])

    def test_sieve_30(self):
        self.assertEqual(sieve(30), [2, 3, 5, 7, 11, 13, 17, 19, 23, 29])

    def test_sieve_2(self):
        self.assertEqual(sieve(2), [2])


class TestSmallPrimeFilter(unittest.TestCase):
    def test_small_primes_pass(self):
        for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]:
            self.assertTrue(passes_small_prime_filter(p), f"{p} should pass")

    def test_composites_filtered(self):
        for n in [4, 6, 8, 9, 10, 14, 21, 25, 27, 49, 77, 91, 121]:
            self.assertFalse(passes_small_prime_filter(n), f"{n} should be filtered")

    def test_large_prime_passes(self):
        self.assertTrue(passes_small_prime_filter(100003))
        self.assertTrue(passes_small_prime_filter(1000003))


class TestMillerRabin(unittest.TestCase):
    def test_small_primes(self):
        for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 97, 101, 103, 107, 109]:
            self.assertTrue(is_probable_prime_mr(p))

    def test_small_composites(self):
        for n in [0, 1, 4, 6, 8, 9, 10, 12, 15, 21, 25, 27, 33]:
            self.assertFalse(is_probable_prime_mr(n))

    def test_large_primes(self):
        primes = [100003, 100019, 100043, 100049, 1000003, 1000033, 10000019, 10000079]
        for p in primes:
            self.assertTrue(is_probable_prime_mr(p), f"{p} should be prime")

    def test_large_composite(self):
        self.assertFalse(is_probable_prime_mr(100003 * 100019))


class TestIsPrimeFast(unittest.TestCase):
    def test_known_primes(self):
        for p in [2, 3, 97, 101, 100003, 1000003, 10000019]:
            self.assertTrue(is_prime_fast(p))

    def test_known_composites(self):
        for n in [0, 1, 4, 6, 100, 1001, 100003 * 100019]:
            self.assertFalse(is_prime_fast(n))


class TestMDigitPrimesFast(unittest.TestCase):
    def test_1_digit(self):
        self.assertEqual(m_digit_primes_fast(1), [2, 3, 5, 7])

    def test_2_digit(self):
        res = m_digit_primes_fast(2)
        self.assertEqual(len(res), 21)
        self.assertEqual(res[0], 11)
        self.assertEqual(res[-1], 97)

    def test_3_digit(self):
        res = m_digit_primes_fast(3)
        self.assertEqual(len(res), 143)
        self.assertEqual(res[0], 101)
        self.assertEqual(res[-1], 997)

    def test_4_digit(self):
        res = m_digit_primes_fast(4)
        self.assertEqual(len(res), 1061)
        self.assertEqual(res[0], 1009)
        self.assertEqual(res[-1], 9973)


class TestMDigitPrimes(unittest.TestCase):
    """m_digit_primes с method='fast' и method='sieve'."""

    def test_1_digit_fast(self):
        self.assertEqual(m_digit_primes(1, method='fast'), [2, 3, 5, 7])

    def test_1_digit_sieve(self):
        self.assertEqual(m_digit_primes(1, method='sieve'), [2, 3, 5, 7])

    def test_2_digit_both_match(self):
        fast = m_digit_primes(2, method='fast')
        sieve_res = m_digit_primes(2, method='sieve')
        self.assertEqual(fast, sieve_res)


class TestPrimeStats(unittest.TestCase):
    def test_stats_tracked(self):
        stats = PrimeStats()
        res = m_digit_primes_fast(2, stats=stats)
        self.assertEqual(len(res), 21)
        self.assertEqual(stats.primes_found, 21)
        # wheel30 yields 18 candidates (24 minus 6 малых простых 11,13,17,19,23,29)
        self.assertEqual(stats.total_candidates, 18)
        self.assertEqual(stats.filtered_by_small_primes, 0)
        self.assertEqual(stats.miller_rabin_calls, 18)


class TestFactorizeByDigits(unittest.TestCase):
    """factorize_by_digits с method='fast' (default)."""
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

    def test_dfs_mode_stack(self):
        self.assertEqual(factorize_by_digits(143, 2, dfs_mode='stack'), [(11, 13)])

    def test_dfs_mode_recursive(self):
        self.assertEqual(factorize_by_digits(143, 2, dfs_mode='recursive'), [(11, 13)])

    def test_dfs_modes_match(self):
        for N, m in [(143, 2), (323, 2), (899, 2), (1147, 2), (1763, 2), (121, 2)]:
            r1 = factorize_by_digits(N, m, dfs_mode='recursive')
            r2 = factorize_by_digits(N, m, dfs_mode='stack')
            self.assertEqual(r1, r2, f"Mismatch for N={N}, m={m}")

    def test_dfs_mode_invalid(self):
        with self.assertRaises(ValueError):
            factorize_by_digits(143, 2, dfs_mode='invalid')


if __name__ == "__main__":
    unittest.main()
