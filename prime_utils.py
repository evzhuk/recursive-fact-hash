import warnings
from dataclasses import dataclass, field

SMALL_PRIMES = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
    31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
]


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


# ---------------------------------------------------------------------------
# Колесо mod 210 = 2×3×5×7: 48 кандидатов из каждых 210
# вместо 8 из 30 — на 37.5% меньше кандидатов чем wheel30
# ---------------------------------------------------------------------------
_WHEEL210_OFFSETS = (
    1, 11, 13, 17, 19, 23, 29, 31, 37, 41,
    43, 47, 53, 59, 61, 67, 71, 73, 79, 83,
    89, 97, 101, 103, 107, 109, 113, 121, 127, 131,
    137, 139, 143, 149, 151, 157, 163, 167, 169, 173,
    179, 181, 187, 191, 193, 197, 199, 209,
)

_WHEEL210_GAPS = tuple(
    _WHEEL210_OFFSETS[(i + 1) % 48] - _WHEEL210_OFFSETS[i]
    if i < 47
    else 210 - _WHEEL210_OFFSETS[47] + _WHEEL210_OFFSETS[0]
    for i in range(48)
)


def _wheel210_start(low: int) -> tuple[int, int]:
    base = (low // 210) * 210
    for i, off in enumerate(_WHEEL210_OFFSETS):
        if base + off >= low:
            return base + off, i
    return base + 210 + _WHEEL210_OFFSETS[0], 0


def _wheel210_candidates(low: int, high: int):
    n, i = _wheel210_start(low)
    gaps = _WHEEL210_GAPS
    while n <= high:
        yield n
        n += gaps[i]
        i = (i + 1) % 48


def m_digit_primes_fast(m: int, stats: PrimeStats | None = None) -> list[int]:
    if m >= 9:
        warnings.warn(
            f"m={m}: генерация простых займёт десятки минут. "
            "Для m>=9 рекомендуется внешняя база или сегментированное решето.",
            RuntimeWarning, stacklevel=2,
        )
    low = 10 ** (m - 1)
    high = 10**m - 1
    result = []

    for p in SMALL_PRIMES:
        if low <= p <= high:
            if stats is not None:
                stats.primes_found += 1
            result.append(p)

    last_small = SMALL_PRIMES[-1]
    for n in _wheel210_candidates(low, high):
        if n <= last_small:
            continue
        if stats is not None:
            stats.total_candidates += 1
            stats.miller_rabin_calls += 1
        if is_probable_prime_mr(n):
            if stats is not None:
                stats.primes_found += 1
            result.append(n)

    return sorted(result)
