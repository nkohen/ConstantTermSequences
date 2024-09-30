def get_coefficient(poly, exponent): # Returns the corresponding coefficient of poly
    R.<t> = LaurentPolynomialRing(ZZ, 1)
    return (poly*t^(-exponent)).constant_coefficient()

def compute_triangle(P, p, num_rows): # Returns a num_rows x num_columns array and the column index of constant terms so that triangle[i][offset + j] = coeff_j(P^i)
    R.<t> = LaurentPolynomialRing(GF(p), 1)

    deg_right = P.degree()
    deg_left = (R(P)(t^-1)).degree()
    num_columns = (num_rows - 1)*(deg_left + deg_right) + 1
    offset = (num_rows - 1)*deg_left # index of constant term
    triangle = []
    poly = 1
    for i in range(num_rows):
        row = [0 for j in range(num_columns)]
        for j in range(-i*deg_left, i*deg_right + 1):
            row[j + offset] = get_coefficient(poly, j)
        triangle.append(row)
        poly = poly*P
    
    return (triangle, offset)

def Lambda(P, p): # If Q(t) is P(t) with all terms whose exponents are not multiples of p deleted, this returns Q(t^{1/p})
    R.<t> = LaurentPolynomialRing(GF(p), 1)
    exps = R(P).exponents()
    d = R(P).dict()
    simp = R(0) # Begin with 0
    for j in range(len(d)):
        exp = exps[j]
        if (exp[0] % p == 0): # Only include terms whose exponent is a multiple of p
            simp += d[exp]*t^(exp[0]/p) # Divide exponent by p
    return simp