"""
Факторизация N = p*q перебором по разрядам справа налево
с использованием дерева вариантов, хэш-структур и рекурсии (или итеративного стека).

Режимы получения простых чисел:
  - 'sieve' (legacy): решето Эратосфена (полный bytearray)
  - 'fast' (default): presieve по первым 10 простым + Miller–Rabin
"""

from array import array
from collections import deque
import time

from prime_utils import is_prime_fast, is_probable_prime_mr, passes_small_prime_filter
from prime_utils import m_digit_primes_fast, PrimeStats


# ---------------------------------------------------------------------------
# Legacy: решето Эратосфена (полная генерация, память ~10^m байт)
# ---------------------------------------------------------------------------

def sieve(limit):
    """Решето Эратосфена: все простые <= limit."""
    if limit > 10**7:
        raise ValueError(f"sieve: limit={limit} > 10^7 требует >100MB — используйте method='fast'")
    is_prime = bytearray(b'\x01') * (limit + 1)
    is_prime[0:2] = b'\x00\x00'
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            step = i
            start = i * i
            is_prime[start:limit+1:step] = b'\x00' * ((limit - start) // step + 1)
    return [i for i, v in enumerate(is_prime) if v]


def m_digit_primes(m, method='fast'):
    """Все простые числа длины ровно m десятичных цифр.

    method='fast' (default) — Miller–Rabin + presieve.
    method='sieve' — решето Эратосфена (legacy).
    """
    if method == 'sieve':
        low = 10 ** (m - 1)
        high = 10**m - 1
        all_primes = sieve(high)
        return [p for p in all_primes if p >= low]
    return m_digit_primes_fast(m)


def build_extend_map(primes, m, timeout_sec=None):
    """
    ext_bits[t, tail] — uint16 битовая маска цифр d (0..9),
    которые можно дописать слева к хвосту tail длины t.
    Бит d установлен если существует простое с таким хвостом и такой следующей цифрой.

    Индексация: t от 1 до m-1, tail от 0 до 10^t - 1.
    Форма массива: (m, 10^(m-1)) — строки по t, столбцы по tail.

    Хранится как плоский array('H') с row_stride = 10^(m-1).
    """
    deadline = time.monotonic() + timeout_sec if timeout_sec else None
    powers = [10**t for t in range(m)]
    cols = powers[m - 1]
    ext_bits = array('H', [0]) * (m * cols)
    for p in primes:
        if deadline is not None and time.monotonic() > deadline:
            raise TimeoutError(f"build_extend_map превысил {timeout_sec}s")
        for t in range(1, m):
            tail = p % powers[t]
            d = (p // powers[t]) % 10
            idx = t * cols + tail
            ext_bits[idx] |= 1 << d
    return ext_bits


def last_digit_pairs(N):
    """
    Пары (a,b) ∈ {1,3,7,9}^2, такие что (a·b) ≡ N (mod 10).
    Простые > 5 могут оканчиваться только на 1,3,7,9 — поэтому
    перебираем только эти цифры.
    """
    r = N % 10
    digits = [1, 3, 7, 9]
    return [(a, b) for a in digits for b in digits if (a * b) % 10 == r]


def factorize_by_digits(N, m, method='fast', dfs_mode='recursive'):
    """
    Вход: N = p·q, m — число десятичных цифр в p и q.
    method='fast' (default) — Miller–Rabin + presieve.
    method='sieve' — решето Эратосфена (legacy).
    dfs_mode='recursive' (default) — рекурсивный DFS.
    dfs_mode='stack' — итеративный DFS через deque.
    Возвращает список пар (p, q) с p ≤ q, p·q = N.
    """
    if not isinstance(N, int) or N <= 1:
        raise ValueError(f"N должно быть целым > 1, получено {N!r}")
    if not isinstance(m, int) or m < 2:
        raise ValueError(f"m должно быть целым >= 2, получено {m!r}")
    if N % 10 == 0:
        raise ValueError(f"N={N} кратно 10 — множители не могут оканчиваться на 0")
    primes = m_digit_primes(m, method=method)
    ext_bits = build_extend_map(primes, m)
    cols = 10**(m - 1)

    # Отбираем стартовые пары последних цифр (единицы),
    # которые встречаются среди m-значных простых
    last_digits_ok = {p % 10 for p in primes}
    start_pairs = [
        (a, b) for a, b in last_digit_pairs(N)
        if a in last_digits_ok and b in last_digits_ok
    ]

    results = set()

    powers = [10**t for t in range(m + 1)]

    if dfs_mode == 'recursive':
        def dfs(t, p_t, q_t):
            if t == m:
                if p_t * q_t == N:
                    results.add((p_t, q_t) if p_t <= q_t else (q_t, p_t))
                return
            mask_p = ext_bits[t * cols + p_t]
            mask_q = ext_bits[t * cols + q_t]
            if mask_p == 0 or mask_q == 0:
                return
            power_t = powers[t]
            mod = powers[t + 1]
            target = N % mod
            for d_p in range(10):
                if not (mask_p >> d_p & 1):
                    continue
                p_next = d_p * power_t + p_t
                for d_q in range(10):
                    if not (mask_q >> d_q & 1):
                        continue
                    q_next = d_q * power_t + q_t
                    if (p_next * q_next) % mod == target:
                        dfs(t + 1, p_next, q_next)
        for a, b in start_pairs:
            dfs(1, a, b)

    elif dfs_mode == 'stack':
        stack = deque((1, a, b) for a, b in start_pairs)
        while stack:
            t, p_t, q_t = stack.pop()

            if t == m:
                if p_t * q_t == N:
                    results.add((p_t, q_t) if p_t <= q_t else (q_t, p_t))
                continue

            mask_p = ext_bits[t * cols + p_t]
            mask_q = ext_bits[t * cols + q_t]
            if mask_p == 0 or mask_q == 0:
                continue

            power_t = powers[t]
            mod = powers[t + 1]
            target = N % mod

            for d_p in range(10):
                if not (mask_p >> d_p & 1):
                    continue
                p_next = d_p * power_t + p_t
                for d_q in range(10):
                    if not (mask_q >> d_q & 1):
                        continue
                    q_next = d_q * power_t + q_t
                    if (p_next * q_next) % mod == target:
                        stack.append((t + 1, p_next, q_next))

    else:
        raise ValueError(f"dfs_mode должен быть 'recursive' или 'stack', получено {dfs_mode!r}")

    return sorted(results)


# ----------------------------------------------------------------------
if __name__ == '__main__':
    # Демонстрация на N = p·q, где p,q — двузначные простые
    test_cases = [
        143,    # 11*13
        323,    # 17*19
        899,    # 29*31
        1147,   # 31*37
        1763,   # 41*43
    ]

    for N in test_cases:
        print(f"\nN = {N}, ищем двузначные множители...")
        res = factorize_by_digits(N, 2)
        if res:
            for p, q in res:
                print(f"  Найдено: {p} * {q} = {N}  (проверка: {p*q == N})")
        else:
            print("  Разложение не найдено.")
