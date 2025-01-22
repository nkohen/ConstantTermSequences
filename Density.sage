def compute_densities(P, Q, p): # Computes the density of each output of ct[P^nQ] mod p
    # Begin by computing time spent on each DFA state, and then sum over states with equal outputs
    # For a reference, see Section 8.2 of Automatic Sequences by Allouche and Shallit
    from DFA import PolyAuto
    R.<t> = LaurentPolynomialRing(GF(p), 1)
    (states, transitions, output_func) = PolyAuto(P, Q, p, 10000)

    mat = [[0 for _ in range(len(states))] for _ in range(len(states))]

    for i in range(len(states)):
        for j in range(p):
            (_, next_state_index) = transitions[i][j]
            mat[i][next_state_index] += 1

    count_mat = Matrix(mat)

    state_values = list(map(output_func, states))

    state_to_value = [[0 for _ in range(p)] for _ in range(len(states))]
    for i in range(len(states)):
        for j in range(p):
            if (state_values[i] == j):
                state_to_value[i][j] = 1

    def kill(x):
        if(x == p):
            return 1
        else:
            return 0
        
    def re(z):
        return z.real()

    (J, P) = Matrix(QQbar, count_mat).jordan_form(transformation=True)

    if (not P.is_invertible()): # For some reason this check is necessary to avoid errors sometimes
        raise Exception(f"Non-invertible matrix unexpected: {P}")

    state_densities = (Matrix([[1] + [0]*(len(states)-1)])*P*J.apply_map(kill)*P.inverse()).apply_map(re)
    return state_densities*Matrix(state_to_value)

def motzkin_zero_density_mod(p): # Returns motzkin zero density according to Corollary 12
    from Sequences import Central_Trinomial, Motzkin_mod

    if(p == 2): # 3-case form only holds for p>2
        return 1/3
    for i in range(p):
        if(Central_Trinomial(i) % p == 0):
            return 1
    
    F = GF(p)
    sgn = F(-3)^((p-1)/2)
    motzkin_zero_count = 0
    likely_zero_count = 0
    unlikely_zero_count = 0

    for i in range(p-1):
        if (F(Central_Trinomial(i)) == sgn*F(Central_Trinomial(i+1))):
            likely_zero_count += 1
        if (F(Central_Trinomial(i)) == F(Central_Trinomial(i+1))):
            unlikely_zero_count += 1

    for i in range(p-2):
        if (Motzkin_mod(i, p)== 0):
            motzkin_zero_count += 1
    
    density = motzkin_zero_count/p + 2*likely_zero_count/((p-1)*(p+1)) + 2*unlikely_zero_count/((p-1)*p*(p+1))
    
    return density

def general_motzkin_zero_density_mod(a, b, p): # Returns motzkin zero density according to Proposition 11
    from Sequences import General_Central_Trinomial, General_Motzkin_mod

    if (p == 2): # 3-case form only holds for p>2
        R.<t> = LaurentPolynomialRing(GF(p), 1)
        return compute_densities(a*t^-1 + b + a*t, 1-t^2, p)[0][0]
    
    if (a % p == 0): # If p | a, then M^{a,b}_n = b^n mod p
        if (b % p == 0):
            return 1
        else:
            return 0

    for i in range(p):
        if(General_Central_Trinomial(a, b, i) % p == 0):
            return 1
    
    F = GF(p)
    sgn = F(b^2 - 4*a^2)^((p-1)/2)
    motzkin_zero_count = 0
    likely_zero_count = 0
    unlikely_zero_count = 0

    for i in range(p-1):
        if (F(b)*F(General_Central_Trinomial(a, b, i)) == sgn*F(General_Central_Trinomial(a, b, i+1))):
            likely_zero_count += 1
        if (F(b)*F(General_Central_Trinomial(a, b, i)) == F(General_Central_Trinomial(a, b, i+1))):
            unlikely_zero_count += 1

    for i in range(p-2):
        if (General_Motzkin_mod(a, b, i, p) == 0):
            motzkin_zero_count += 1
    
    density = motzkin_zero_count/p + 2*likely_zero_count/((p-1)*(p+1)) + 2*unlikely_zero_count/((p-1)*p*(p+1))
    
    return density

def generic_linear_zero_density_mod(Q, a, b, c, d, p): # Returns motzkin zero density according to the process of Section 3.1
    from Sequences import Central_Trinomial, Constant_Term_mod

    F = GF(p)
    R.<t> = LaurentPolynomialRing(F, 1)

    if(p == 2): # 3-case form only holds for p>2
        return compute_densities(t^-1 + 1 + t, Q, p)[0][0]
    for i in range(p):
        if(Central_Trinomial(i) % p == 0):
            return 1

    sgn = F(-3)^((p-1)/2)
    zero_count = 0
    likely_zero_count = 0
    unlikely_zero_count = 0
    zs = []
    ls = []
    us = []

    for i in range(p-1):
        if (Constant_Term_mod(t^-1 + 1 + t, Q, i, p) == 0):
            zs.append(i)
            zero_count += 1
        if (a*F(Central_Trinomial(i)) == b*sgn*F(Central_Trinomial(i+1))):
            ls.append(i)
            likely_zero_count += 1
        if (c*F(Central_Trinomial(i)) == d*F(Central_Trinomial(i+1))):
            us.append(i)
            unlikely_zero_count += 1
    
    density = zero_count/p + likely_zero_count/((p-1)*(p+1)) + unlikely_zero_count/((p-1)*p*(p+1))
    return density