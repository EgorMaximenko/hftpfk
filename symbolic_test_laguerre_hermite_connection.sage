
def fourier_laguerre_of_square_sum(n, a, u, xi):
    expr0 = exp(- I * u * xi) * exp(- u ** 2 / 2) * laguerre(n, u ** 2 + a ** 2)
    expr1 = integral(expr0, u, -Infinity, Infinity)
    expr2 = expr1 / sqrt(2 * pi)
    return expr2.full_simplify()


def pol_fourier_laguerre_of_square_sum(n, a, u, xi):
    expr0 = fourier_laguerre_of_square_sum(n, a, u, xi)
    factor0 = exp(xi ** 2 / 2)
    factor1 = factorial(n)
    expr1 = factor0 * factor1 * expr0
    return expr1.full_simplify()


def product_hermite(n, a, u, xi):
    expr0 = hermite(n, (xi + a) / sqrt(2))
    expr1 = hermite(n, (xi - a) / sqrt(2))
    expr2 = exp(- xi ** 2 / 2) / ((2 ** n) * factorial(n))
    expr3 = expr0 * expr1 * expr2
    return expr3.full_simplify()


def pol_product_hermite(n, a, u, xi):
    expr0 = hermite(n, (xi + a) / sqrt(2))
    expr1 = hermite(n, (xi - a) / sqrt(2))
    expr2 = 1 / (2 ** n)
    expr3 = expr0 * expr1 * expr2
    return expr3.full_simplify()


def fourier_product_hermite(n, a, u, xi):
    factor0 = product_hermite(n, a, u, xi)
    factor1 = exp(I * u * xi)
    expr0 = factor0 * factor1
    expr1 = integral(expr0, xi, -Infinity, Infinity)
    expr2 = expr1 / sqrt(2 * pi)
    return expr2.full_simplify()


def laguerre_of_square_sum(n, a, u, xi):
    expr0 = exp(- u ** 2 / 2) * laguerre(n, u ** 2 + a ** 2)
    return expr0.full_simplify()


def big_test_fourier_laguerre_of_square_sum(nmax):
    a, u, xi = var(['a', 'u', 'xi'])
    assume(a > 0)
    big_result = True
    for n in range(nmax + 1):
        result0 = fourier_laguerre_of_square_sum(n, a, u, xi)
        result1 = product_hermite(n, a, u, xi)
        result2 = result0 - result1
        result3 = result2.full_simplify()
        result4 = result3.is_zero()
        big_result = big_result and result4
        print('n = ', n)
        print('  lhs = ', result0)
        print('  rhs = ', result1)
        print('  lhs == rhs?', result4, '\n')
    return big_result


def big_test_fourier_product_hermite(nmax):
    a, u, xi = var(['a', 'u', 'xi'])
    assume(a > 0)
    big_result = True
    for n in range(nmax + 1):
        result0 = laguerre_of_square_sum(n, a, u, xi)
        result1 = fourier_product_hermite(n, a, u, xi)
        result2 = result0 - result1
        result3 = result2.full_simplify()
        result4 = result3.is_zero()
        big_result = big_result and result4
        print('n = ', n)
        print('  lhs = ', result0)
        print('  rhs = ', result1)
        print('  lhs == rhs?', result4, '\n')
    return big_result
    

def big_test_pol_fourier_laguerre_of_square_sum(nmax):
    a, u, xi = var(['a', 'u', 'xi'])
    assume(a > 0)
    big_result = True
    for n in range(nmax + 1):
        result0 = pol_fourier_laguerre_of_square_sum(n, a, u, xi)
        result1 = pol_product_hermite(n, a, u, xi)
        result2 = result0 - result1
        result3 = result2.full_simplify()
        result4 = result3.is_zero()
        big_result = big_result and result4
        print('n = ', n)
        print('  lhs = ', result0)
        print('  rhs = ', result1)
        print('  lhs == rhs?', result4, '\n')
    return big_result


#print(big_test_fourier_laguerre_of_square_sum(10))
print(big_test_fourier_product_hermite(10))
#print(big_test_pol_fourier_laguerre_of_square_sum(10))

