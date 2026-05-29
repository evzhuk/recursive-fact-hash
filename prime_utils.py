from dataclasses import dataclass, field

SMALL_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]

MR_BASES = [2, 325, 9375, 28178, 450775, 9780504, 1795265022]


@dataclass
class PrimeStats:
    total_candidates: int = 0
    filtered_by_small_primes: int = 0
    miller_rabin_calls: int = 0
    primes_found: int = 0


def passes_small_prime_filter(n: int) -> bool:
    if n < 2:
        return False
    for p in SMALL_PRIMES:
        if n == p:
            return True
        if n % p == 0:
            return False
    return True


def _mr_bases_for_n(n: int) -> list[int]:
    if n < 2047:
        return [2]
    if n < 1373653:
        return [2, 3]
    if n < 9080191:
        return [31, 73]
    if n < 25326001:
        return [2, 3, 5]
    if n < 3215031751:
        return [2, 3, 5, 7]
    return [2, 325, 9375, 28178, 450775, 9780504, 1795265022]


def is_probable_prime_mr(n: int) -> bool:
    if n < 2:
        return False
    for p in SMALL_PRIMES:
        if n == p:
            return True
        if n % p == 0:
            return False

    d = n - 1
    s = 0
    while d & 1 == 0:
        d >>= 1
        s += 1

    for a in _mr_bases_for_n(n):
        if a >= n:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True


def is_prime_fast(n: int) -> bool:
    if n < 2:
        return False
    for p in SMALL_PRIMES:
        if n == p:
            return True
        if n % p == 0:
            return False
    return is_probable_prime_mr(n)


def m_digit_primes_fast(m: int, stats: PrimeStats | None = None) -> list[int]:
    low = 10 ** (m - 1)
    high = 10**m - 1
    result = []
    for n in range(low, high + 1):
        if stats is not None:
            stats.total_candidates += 1
        if not passes_small_prime_filter(n):
            if stats is not None:
                stats.filtered_by_small_primes += 1
            continue
        if stats is not None:
            stats.miller_rabin_calls += 1
        if is_probable_prime_mr(n):
            if stats is not None:
                stats.primes_found += 1
            result.append(n)
    return result
