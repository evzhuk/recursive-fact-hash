import time
import sys
from factorize_by_digits import m_digit_primes, factorize_by_digits, build_extend_map
from prime_utils import PrimeStats, m_digit_primes_fast


def bench(label, fn):
    start = time.perf_counter()
    result = fn()
    elapsed = time.perf_counter() - start
    print(f"  {label:<30s} {elapsed:>10.4f}s")
    return result, elapsed


def run_fast(m, N):
    stats = PrimeStats()

    primes, t1 = bench("  m_digit_primes_fast", lambda: m_digit_primes_fast(m, stats=stats))
    ext, t2 = bench("  build_extend_map", lambda: build_extend_map(primes, m))

    t3s = time.perf_counter()
    res = _factorize_with_precomputed(primes, ext, N, m)
    t3 = time.perf_counter() - t3s
    print(f"  {'DFS only':<30s} {t3:>10.4f}s")
    print(f"  {'→ result':<30s} {str(res):>10}")
    print(f"  {'→ total precompute+DFS':<30s} {t1+t2+t3:>10.4f}s")

    print(f"\n  Statistics:")
    print(f"    total_candidates        {stats.total_candidates:>10}")
    print(f"    filtered_by_small_primes {stats.filtered_by_small_primes:>10}")
    print(f"    miller_rabin_calls      {stats.miller_rabin_calls:>10}")
    print(f"    primes_found            {stats.primes_found:>10}")
    if stats.total_candidates > 0:
        pct = 100 * stats.filtered_by_small_primes / stats.total_candidates
        print(f"    small-prime filter rate  {pct:>9.1f}%")

    return primes, res


def run_sieve(m, N):
    primes, t1 = bench("  m_digit_primes (sieve)", lambda: m_digit_primes(m, method='sieve'))
    ext, t2 = bench("  build_extend_map", lambda: build_extend_map(primes, m))

    t3s = time.perf_counter()
    res = _factorize_with_precomputed(primes, ext, N, m)
    t3 = time.perf_counter() - t3s
    print(f"  {'DFS only':<30s} {t3:>10.4f}s")
    print(f"  {'→ result':<30s} {str(res):>10}")
    print(f"  {'→ total precompute+DFS':<30s} {t1+t2+t3:>10.4f}s")

    return primes, res


def _factorize_with_precomputed(primes, ext, N, m):
    last_digits_ok = {p % 10 for p in primes}
    start_pairs = []
    r = N % 10
    digits = [1, 3, 7, 9]
    for a in digits:
        if a not in last_digits_ok:
            continue
        for b in digits:
            if b not in last_digits_ok:
                continue
            if (a * b) % 10 == r:
                start_pairs.append((a, b))

    results = set()

    def dfs(t, p_t, q_t):
        if t == m:
            if p_t * q_t == N:
                results.add((p_t, q_t) if p_t <= q_t else (q_t, p_t))
            return

        p_candidates = ext.get((t, p_t))
        q_candidates = ext.get((t, q_t))
        if p_candidates is None or q_candidates is None:
            return

        mod = 10 ** (t + 1)
        target = N % mod

        for d_p in p_candidates:
            p_next = d_p * (10**t) + p_t
            for d_q in q_candidates:
                q_next = d_q * (10**t) + q_t
                if (p_next * q_next) % mod == target:
                    dfs(t + 1, p_next, q_next)

    for a, b in start_pairs:
        dfs(1, a, b)

    return sorted(results)


def benchmark(m):
    N = known_N[m]

    print(f"\n{'='*60}")
    print(f"  m = {m}  ({m}-digit factors, N ≈ {N})")
    print(f"{'='*60}")

    print(f"\n  --- Miller–Rabin (fast) ---")
    primes_fast, res_fast = run_fast(m, N)

    if m <= 7:
        print(f"\n  --- Sieve of Eratosthenes (legacy) ---")
        primes_sieve, res_sieve = run_sieve(m, N)
        match = (primes_fast == primes_sieve) and (res_fast == res_sieve)
        print(f"\n  Results match: {match}")
    else:
        print(f"\n  --- Sieve (legacy) ---")
        print(f"  Пропущен: m={m} требует >100MB для sieve bytearray")


known_N = {
    2: 143,
    3: 179 * 191,
    4: 1009 * 1013,
    5: 10007 * 10009,
    6: 100003 * 100019,
    7: 1000003 * 1000033,
    8: 10000019 * 10000079,
}

if __name__ == "__main__":
    max_m = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    for m in range(2, max_m + 1):
        if m not in known_N:
            print(f"\nm={m}: пропускаем (нет тестового N)")
            continue
        try:
            benchmark(m)
        except (MemoryError, KeyboardInterrupt) as e:
            print(f"  ОШИБКА: {e}")
            break
