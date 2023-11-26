import math
import random
import sympy


# 1a) Euler function (find quantity of number which are less and mutually prime with current number)
def Euler_function(n):
    print(f"Euler function for {n}:")
    mult = 1
    i = 2
    while i * i <= n:
        p = 1
        while n % i == 0:
            n = n / i
            p = p * i
        p = p / i
        if p >= 1:
            mult = mult * p * (i - 1)
        i = i + 1
    n = n - 1

    if n == 0:
        print(int(mult))
    else:
        print(int(n * mult))


# 1b) Möbius function
def Mobius_function(n):
    print(f"Möbius function for {n}:")

    if n == 1:
        print(1)
    p = 0
    for i in range(1, n + 1):
        if (n % i == 0 and
                sympy.isprime(i)):
            if n % (i * i) == 0:
                print(0)
            else:
                p = p + 1
    if p % 2 != 0:
        print(-1)
    else:
        print(1)


# 2) Chinese remainder theorem
def Chinese_remainder_theorem(num, rem):
    print(f"Chinese remainder theorem for system:")
    for i in range(len(rem)):
        print(f"N = {num[i]} (mod {rem[i]}) ")

    M = 1
    for m in rem:
        M *= m

    s = [0] * len(rem)
    x = 0
    for i in range(len(rem)):
        y = M / rem[i] % rem[i]
        j = 1
        while j < rem[i] and s[i] == 0:
            if math.gcd(j, rem[i]) == 1:
                if y * j % rem[i] == 1:
                    s[i] = j
            j += 1
        x += num[i] * (M / rem[i]) * s[i]
    result = int(x % M)
    print("N = ", result)


# Additional function for calculating Legendre Symbols
def part_of_div(q, p):
    div = 1
    while not sympy.isprime(q) and q != 1:
        prime_multiplier = int(factorization(q))
        div *= part_of_div(prime_multiplier, p)
        q = int(q / prime_multiplier)
        div *= part_of_div(q, p)
        return div
    else:
        if q == 2:
            if p % 8 == 1 or p % 8 == 7:
                return div
            elif p % 8 == 3 or p % 8 == 5:
                return -div
        elif int(math.sqrt(q)) ** 2 == q:
            return div
        else:
            if p % 4 == 3 and q % 4 == 3:
                p = p % q
                return -div * part_of_div(p, q)
            else:
                p = p % q
                return div * part_of_div(p, q)


# 3a) Legendre Symbols calculating (p should be prime, positive and odd)
def Legendre_symbols(a, p):
    a = a % p
    if a == 0:
        div = 0
    else:
        div = 1
        div *= part_of_div(a, p)
    return div


# 3a) Legendre Symbols print
def Legendre_symbols_print(a, p):
    print(f"Legendre symbols for {a}/{p}: ")
    if not sympy.isprime(p) or p == 2 or p <= 0:
        print("Legendre symbol is undefined")
        return 0

    div = Legendre_symbols(a, p)
    print(div)
    return 0


# 4b) Jacobi Symbols (p should be positive and odd)
def Jacobi_symbols(a, p):
    print(f"Jacobi symbols for {a}/{p}: ")
    if p == 2:
        print("Jacobi symbol is undefined")
        return 0

    a = a % p
    if a == 0:
        div = 0
    else:
        div = 1
        while not sympy.isprime(p):
            prime_multiplier = int(factorization(p))
            div *= Legendre_symbols(a, prime_multiplier)
            p = int(p / prime_multiplier)
        div *= Legendre_symbols(a, p)
    print(div)
    return 0


# Additional function for factorization (Polynomial)
def g(x, n, c):
    return (x ** 2 + c) % n


# Pollard's rho algorithm for factorization
def factorization(n):
    x = 2
    y = x
    m = 2
    d = 1
    c = 1

    while 1:
        if d == 1:
            if int(math.sqrt(n)) ** 2 == n:
                if not sympy.isprime(int(math.sqrt(n))):
                    d = factorization(int(math.sqrt(n)))
                else:
                    d = math.sqrt(n)
                return d
            x = g(x, n, c)
            y = g(g(y, n, c), n, c)
            d = math.gcd(abs(x - y), n)
        else:
            if d == n:
                m = m + 1
                x = m + 1
                if x > n - 1:
                    return 0
                y = x
                d = 1
            else:
                return d


# 4) Pollard's rho algorithm for factorization print
def factorization_print(n):
    d = factorization(n)
    print(f"Pollard's rho algorithm (factorization) for {n}:")
    if d == 0:
        print("This number is prime")
    else:
        print('{0:.0f}'.format(d))


# 5) Baby step - Giant step algorithm for logarithms
def logarithmization(a, b, n):  # a^x=b
    print(f"Baby step - Giant step (logarithmization) {a}^x(mod {n})={b}:")
    k = round(math.sqrt(n)) + 1
    baby = [0] * k
    giant = [0] * k

    for i in range(k):
        baby[i] = pow(a, i) % n

    for j in range(k):
        giant[j] = (b * (pow(pow(a, k * (n - 2)) % n, j) % n)) % n
        if giant[j] in baby:
            print(baby.index(giant[j]) + j * k % n)
            return 0
    return 0


def mult_complex(k, j, p, w2):
    real = sympy.re(k[j])
    imag = sympy.im(k[j]) / sympy.sqrt(abs(w2))
    if abs(real) >= p:
        real = abs(real) % p
        if sympy.re(k[j]) < 0:
            real = -real
    if abs(imag) >= p:
        imag = abs(imag) % p
        if sympy.im(k[j]) < 0:
            imag = -imag
    k[j] = real + imag * sympy.sqrt(abs(w2)) * sympy.I

# 6) Cipolla's algorithm
def Cipollas_algorithm(n, p):
    print(f"Cipolla's algorithm for x = {n}(mod {p}):")
    if Legendre_symbols(n, p) != 1:
        print(f"We can't find solution")
    a = 2
    w2 = a * a - n
    w2_rem = w2 % p

    while Legendre_symbols(w2_rem, p) != -1:
        a += 1

    power = int((p + 1) / 2)
    k = [0] * (power+1)
    k[1] = a + sympy.sqrt(w2)
    j = 1
    while j * 2 < power:
        j *= 2
        k[j] = sympy.expand(k[int(j/2)]**2)
        mult_complex(k, j, p, w2)
    sub = power - j
    l = j
    while sub > 0:
        if (k[sub] != 0):
            k[j+sub] = sympy.expand(k[j]*k[sub])
            mult_complex(k, j+sub, p, w2)
            sub = power - (j + sub)
        else:
            while l > sub:
                l /= 2
            sub = sub - l
            l = int(l)
            k[j+l] = sympy.expand(k[j]*k[l])

            j += l
            mult_complex(k, j, p, w2)
    pos_res = sympy.re(k[power])
    neg_res = -pos_res
    print(f"Solutions: {pos_res}, {neg_res}")


# 7) Miller-Rabin algorithm (Primality test)
def Miller_Rabin_algorithm(n, k):
    if n % 2 == 0:
        return False
    m = (n - 1) / 2
    t = 1
    while m % 2 == 0:
        m /= 2
        t += 1
    for i in range(1, k):
        a = random.randint(3, n-4)
        u = a ** int(m) % n
        if u != 1:
            j = 1
            while u != n - 1 and j < t:
                u = u ** 2 % n
                j += 1
            if u != n - 1:
                return False
        else:
            return True
    return True


def Miller_Rabin_algorithm_print(n, k):
    print(f"Miller-Rabin algorithm for {n}:")
    result = Miller_Rabin_algorithm(n, k)
    if result:
        print(f"{n} is prime")
    else:
        print(f"{n} is not prime")

# 8) RSA Encryption
def RSA_Encryption(m):
    print(f"RSA-Encryption for {m}:")
    p = random.randint(3, 1000)
    while not Miller_Rabin_algorithm(p, 10):
        p = random.randint(3, 1000)

    q = random.randint(3, 1000)
    while not Miller_Rabin_algorithm(q, 10) or q == p:
        q = random.randint(3, 1000)

    n = p * q
    lam = sympy.lcm(p - 1, q - 1)

    e = random.randint(2, lam)
    while not Miller_Rabin_algorithm(e, 10) or sympy.gcd(e, lam) != 1:
        e = random.randint(2, lam)

    d = random.randint(2, lam)
    while e * d % lam != 1:
        d = random.randint(2, lam)

    c = m ** e % n
    print(c)
    return c, n, d

# 8) RSA Decryption
def RSA_Decryption(c, n, d):
    print(f"RSA-Decryption for {c}:")
    m = c ** d % n
    print(m)


def add_points(P1, P2, p, a):
    x1, y1 = P1
    x2, y2 = P2

    if x1 == x2 and y1 == y2:
        beta = (3 * x1 * x2 + a) * pow(2 * y1, -1, p)
    else:
        beta = (y2 - y1) * pow(x2 - x1, -1, p)

    x3 = (beta * beta - x1 - x2) % p
    y3 = (beta * (x1 - x3) - y1) % p
    return x3, y3


def double_and_add_method(G, k, p, a):
    k_binary = bin(k)[2:]
    point1 = G
    point2 = G

    for i in range(1, len(k_binary)):
        current_bit = k_binary[i: i + 1]
        if current_bit == '1':
            point1 = add_points(point1, point2, p, a)
        else:
            point2 = add_points(point2, point2, p, a)
    return point1


# 9) Elliptic Curve ElGamal Encryption
def El_Gamal(a, p, n, t, P):
    print(f"Elliptic Curve ElGamal Encryption:")
    k = random.randint(1, n-1)
    r = random.randint(1, n-1)

    Y = double_and_add_method(P, k, p, a)
    Alice_point = double_and_add_method(P, t, p, a)
    print("Alice point: ", Alice_point)

    g = double_and_add_method(P, r, p, a)
    h = add_points(Alice_point, double_and_add_method(Y, r, p, a), p, a)

    s = double_and_add_method(g, k, p, a)
    x, y = s
    y = -y
    s = (x, y)

    Bob_point = add_points(s, h, p, a)
    print("Bob point: ", Bob_point)


if __name__ == '__main__':
    Euler_function(36)
    # Mobius_function(18644564656756586767)
    # Chinese_remainder_theorem([6, 13, 28], [25, 36, 49])
    # Legendre_symbols_print(7, 13)
    # Jacobi_symbols(91, 103)
    # factorization_print(121)
    # logarithmization(2, 2, 41)
    # Cipollas_algorithm(7, 13)
    # Miller_Rabin_algorithm_print(36, 10)
    # c, n, d = RSA_Encryption(7643)
    # RSA_Decryption(c, n, d)
    # El_Gamal(6564647656756756567564, 1564557676556756555423764748464879078528375647687567345687676545638475435434737,
    #          115792089237316195423570985008687907852837564279074904382605163141518161494337,
    #          987877,
    #          (55066263022277343669578718895168534326250603453777594175500187360389116729240,
    #           32670510020758816978083085130507043184471273380659243275938904335757337482424)
    #          )

