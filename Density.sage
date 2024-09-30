def compute_densities(P, Q, p): # Computes the density of each output of ct[P^nQ] mod p
    # Begin by computing time spent on each DFA state, and then sum over states with equal outputs
    import DFA
    R.<t> = LaurentPolynomialRing(GF(p), 1)
    (states, transitions, output_func) = DFA.PolyAuto(P, Q, p, 10000)

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
    state_densities = (Matrix([[1] + [0]*(len(states)-1)])*P*J.apply_map(kill)*P^(-1)).apply_map(re)
    return state_densities*Matrix(state_to_value)

# TODO: Verify that this is in agreement with the Corollary and check against above def
def motzkin_zero_density_mod(p): # Returns motzkin zero density according to Corollary 12
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
        if (F(Central_Trinomial(i)) == F(Central_Trinomial(i+1))):
            unlikely_zero_count += 1
        if (sgn*F(Central_Trinomial(i)) == F(Central_Trinomial(i+1))):
            likely_zero_count += 1

    for i in range(p-2):
        if (Motzkin_mod(i, p)== 0):
            motzkin_zero_count += 1
    
    density = motzkin_zero_count/p + 2*likely_zero_count/((p-1)*(p+1)) + 2*unlikely_zero_count/((p-1)*p*(p+1))
    return (density, motzkin_zero_count, likely_zero_count, unlikely_zero_count)
