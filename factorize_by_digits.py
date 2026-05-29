"""
Факторизация N = p*q перебором по разрядам справа налево
с использованием дерева вариантов, хэш-структур и рекурсии.

Режимы получения простых чисел:
  - 'sieve' (legacy): решето Эратосфена (полный bytearray)
  - 'fast' (default): presieve по первым 10 простым + Miller–Rabin
"""

from prime_utils import is_prime_fast, is_probable_prime_mr, passes_small_prime_filter
from prime_utils import m_digit_primes_fast, PrimeStats


# ---------------------------------------------------------------------------
# Legacy: решето Эратосфена (полная генерация, память ~10^m байт)
# ---------------------------------------------------------------------------

def sieve(limit):
    """Решето Эратосфена: все простые <= limit."""
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


def build_extend_map(primes, m):
    """
    extend_map[(t, tail)] -> set[int]  (цифры d, дописываемые слева)

    Для каждого простого P = d_{m-1}...d_1 (d_1 — единицы):
      для t = 1..m-1:
        tail = P % 10^t
        d   = (P // 10^t) % 10   (следующая слева цифра)
        extend_map[(t, tail)].add(d)

    Если tail не может быть хвостом никакого простого,
    то ключа в словаре нет — ветка обрезается.
    """
    # Для каждого простого числа разбиваем его на хвост (младшие t цифр)
    # и цифру d, которая может стоять слева от этого хвоста.
    # ext[(t, tail)] — множество цифр d, при которых d*10^t + tail — начало простого.
    ext = {}
    for p in primes:
        for t in range(1, m):
            tail = p % (10**t)       # младшие t цифр числа p
            d = (p // (10**t)) % 10  # следующая цифра слева от хвоста
            key = (t, tail)
            if key not in ext:
                ext[key] = set()
            ext[key].add(d)
    return ext


def last_digit_pairs(N):
    """
    Пары (a,b) ∈ {1,3,7,9}^2, такие что (a·b) ≡ N (mod 10).
    Простые > 5 могут оканчиваться только на 1,3,7,9 — поэтому
    перебираем только эти цифры.
    """
    r = N % 10
    digits = [1, 3, 7, 9]
    return [(a, b) for a in digits for b in digits if (a * b) % 10 == r]


def factorize_by_digits(N, m, method='fast'):
    """
    Вход: N = p·q, m — число десятичных цифр в p и q.
    method='fast' (default) — Miller–Rabin + presieve.
    method='sieve' — решето Эратосфена (legacy).
    Возвращает список пар (p, q) с p ≤ q, p·q = N.
    """
    primes = m_digit_primes(m, method=method)
    ext = build_extend_map(primes, m)

    # Отбираем стартовые пары последних цифр (единицы),
    # которые встречаются среди m-значных простых
    last_digits_ok = {p % 10 for p in primes}
    start_pairs = [
        (a, b) for a, b in last_digit_pairs(N)
        if a in last_digits_ok and b in last_digits_ok
    ]

    results = set()

    def dfs(t, p_t, q_t):
        """t — сколько разрядов уже зафиксировано,
        p_t, q_t — части чисел, состоящие из младших t цифр."""
        # Все m разрядов заполнены — проверяем точное равенство
        if t == m:
            if p_t * q_t == N:
                results.add((p_t, q_t) if p_t <= q_t else (q_t, p_t))
            return

        # Какие цифры можно дописать слева к p_t и q_t?
        p_candidates = ext.get((t, p_t))
        q_candidates = ext.get((t, q_t))
        if p_candidates is None or q_candidates is None:
            return  # тупик: такой хвост не может быть началом простого

        # Каким должен быть остаток произведения по модулю 10^(t+1)
        mod = 10 ** (t + 1)
        target = N % mod

        # Перебираем все пары цифр-кандидатов, проверяя условие mod
        for d_p in p_candidates:
            p_next = d_p * (10**t) + p_t
            for d_q in q_candidates:
                q_next = d_q * (10**t) + q_t
                if (p_next * q_next) % mod == target:
                    dfs(t + 1, p_next, q_next)

    # Запускаем DFS от каждой допустимой пары последних цифр
    for a, b in start_pairs:
        dfs(1, a, b)

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
