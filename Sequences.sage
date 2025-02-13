from functools import reduce

def Constant_Term(P, Q, n): # Naively computes ct[P^nQ]
    return ((P^n)*Q).constant_coefficient()

def Constant_Term_mod_slow(P, Q, n, p): # Naively computes ct[P^nQ] mod p
    R.<t> = LaurentPolynomialRing(IntegerModRing(p), 1)
    return ((R(P)^n)*R(Q)).constant_coefficient()

def Constant_Term_mod(P, Q, n, p): # Efficiently computes ct[P^nQ] mod p
    if (n == 0):
        return Q.constant_coefficient() % p

    import PolyUtil
    r = n % p
    new_n = (n - r)/p
    new_Q = PolyUtil.Lambda(Q*P^r, p)
    return Constant_Term_mod(P, new_Q, new_n, p)

motzkin_array = [1,1] # Cached Motzkin sequence
def Motzkin(n): # A001006
    current_len = len(motzkin_array)
    if (n >= current_len): # Extend the cache if needed
        for i in range(n - current_len + 1):
            # Compute next value using (n+2)M_n = (2n+1)M_{n-1} + 3(n-1)M_{n-2}
            motzkin_array.append(((2*(i+current_len) + 1)*motzkin_array[len(motzkin_array) - 1] + (3*(i+current_len) - 3)*motzkin_array[len(motzkin_array) - 2])/(i+2+current_len))  
    
    return motzkin_array[n]

central_array = [1,1] # Cached Central Trinomial Sequence
def Central_Trinomial(n): # A002426
    current_len = len(central_array)
    if (n >= current_len): # Extend the cache if needed
        for i in range(n - current_len + 1):
            # Compute next value using nT_n = (2n-1)T_{n-1} + 3(n-1)T_{n-2}
            central_array.append(((2*(i+current_len) - 1)*central_array[len(central_array) - 1] + (3*(i+current_len) - 3)*central_array[len(central_array) - 2])/(i+current_len))  
    
    return central_array[n]

def General_Central_Trinomial(a, b, n):
    if (n == 0):
        return 1

    prev = 1
    next = b
    for i in range(2, n+1):
        # Compute next value using nT_n = b(2n-1)T_{n-1} - (b^2 - 4a^2)(n-1)T_{n-2}
        new_next = (b*(2*i - 1)*next - (b^2 - 4*a^2)*(i - 1)*prev)/i
        prev = next
        next = new_next
    
    return next

def General_Motzkin(a, b, n):
    if (n == 0):
        return 1

    prev = 1
    next = b
    for i in range(2, n+1):
        # Compute next value using (n+2)M_n = b(2n+1)M_{n-1} - (b^2 - 4a^2)(n-1)M_{n-2}
        new_next = (b*(2*i + 1)*next - (b^2 - 4*a^2)*(i - 1)*prev)/(i+2)
        prev = next
        next = new_next
    
    return next

def Central_Trinomial_mod(n, p): # Fast computation for A002426(n) mod p using the Lucas Congruence
    if (n == 0):
        return 1
    
    R = IntegerModRing(p)

    def central_h(i): # Returns A002426 mod p computed without tricks (used for i < p)
        return R(Central_Trinomial(i))
    def mult_p(a, b): # Returns a*b mod p
        return R(a)*R(b)
    
    return reduce(mult_p, list(map(central_h, Integer(n).digits(p)))) # T_n is congruent to the product of T_{n_i} where n_i are the digits of n in base p

def General_Central_Trinomial_mod(a, b, n, p): # Fast computation for T^{a,b}(n) mod p using the Lucas Congruence
    if (n == 0):
        return 1
    
    R = IntegerModRing(p)

    def central_h(i): # Returns A002426 mod p computed without tricks (used for i < p)
        return R(General_Central_Trinomial(a, b, i))
    def mult_p(a, b): # Returns a*b mod p
        return R(a)*R(b)
    
    return reduce(mult_p, list(map(central_h, Integer(n).digits(p)))) # T_n is congruent to the product of T_{n_i} where n_i are the digits of n in base p

def Motzkin_mod(n, p): # Fast computation for A001006(n) mod p using fast computation for A002426
    R = IntegerModRing(p)
    return (R(2)^(-1))*(3*Central_Trinomial_mod(n, p) + 2*Central_Trinomial_mod(n+1, p) - Central_Trinomial_mod(n+2, p)) # 2M_n = 3T_n + 2T_{n+1} - T_{n+2}

def General_Motzkin_mod(a, b, n, p): # Fast computation for M^{a,b} mod p using fast computation for T^{a,b}
    R = IntegerModRing(p)
    if (R(a) == 0):
        return R(b)^n
    return (R(2*a^2)^(-1))*((4*a^2 - b^2)*General_Central_Trinomial_mod(a, b, n, p) + 2*b*General_Central_Trinomial_mod(a, b, n+1, p) - General_Central_Trinomial_mod(a, b, n+2, p)) # 2a^2M_n = (4a^2 - b^2)T_n + 2bT_{n+1} - T_{n+2}

def next_Motzkin_mod(prev, prev_n, p): # Even faster computation for A001006(prev_n + 1) mod p if A001006(prev_n) mod p is known
    prev_n0 = prev_n % p
    # If the least significant digit is < p-4, then we can undo the effect of this digit in M_{prev_n} computed using the Motzkin DFA and apply the effect of the successor digit
    if (prev_n0 < p-4):
        R = IntegerModRing(p)
        diff_mult = R(Motzkin_mod(prev_n0 + 1, p))*R(Motzkin_mod(prev_n0, p))^(-1)
        return prev*diff_mult
    else: # Otherwise, just default to computing directly
        return Motzkin_mod(prev_n + 1, p)

def Motzkin_range_mod(index, size, p): # Returns [M_index mod p, ..., M_{index+size-1} mod p]
    prev = Motzkin_mod(index, p)
    result = [prev]
    for i in range(size-1):
        prev = next_Motzkin_mod(prev, index + i, p)
        result.append(prev)
    
    return result

def Motzkin_mod_p_mat(n, p): # Computes A001006(n) mod p using the regular representation
    R.<t> = LaurentPolynomialRing(ZZ, 1)
    f = t + 1 + t^(-1)

    def get_mat(k): # Computes the matrix-valued morphism in the regular representation of M_n mod p
        corner = 0
        if k == p-1:
            corner = 1
        poly = f^k
        center = poly.constant_coefficient() % p
        side = 0
        if k > 0:
            side = poly.coefficients()[k+1] % p
        return Matrix([[corner, 0, 0], [side, center, side], [0, 0, corner]])

    def get_seed(k): # Computes the column vector in the regular representation of M_n mod p
        r = k % p

        mr = 1
        if (r > 1):
            poly = f^r
            mr = (poly.coefficients()[r] - poly.coefficients()[r+2]) % p

        if r < p-2:
            return Matrix([[0], [mr], [0]])
        elif r == p-2:
            return Matrix([[0], [mr], [-1]])
        else:
            return Matrix([[0], [mr], [1]])
    
    # Compute using the regular representation of M_n mod p
    state = get_seed(n)
    for digit in Integer(n//p).digits(p):
        state = get_mat(digit)*state
    return state[1][0] # Equivalent to multiplying by the row vector [1, 0, ..., 0]