from PolyUtil import get_coefficient

'''
In this file, when polynomials are interpreted as row vectors they are of
dimension 2*max_deg + 1, and their entries run from the coefficient
of t^(-max_deg) to that of t^(max_deg).
'''

def compute_mat_for_poly(poly, max_deg, p): # Computes the matrix for Q -> Lambda_p(poly*Q)
    return Matrix(GF(p), 2*max_deg + 1,
        lambda i, j: get_coefficient(poly, (max_deg - i)-(max_deg - j)*p)
    )

def compute_mat(P, k, max_deg, p): # Computes the matrix for Q -> Lambda_p((P^k)*Q)
    return compute_mat_for_poly(P^k, max_deg, p)

def poly_to_vec(poly, max_deg): # Returns poly as a vector
    # A more efficient version of Matrix(1, 2*max_deg + 1, lambda _, i : get_coefficient(poly, i-max_deg))
    exps = poly.exponents()
    d = poly.dict()

    current_exp_ind = 0
    vec = []
    for i in range(max_deg, -max_deg-1, -1):
        if (current_exp_ind == len(exps)):
            vec.append(0)
        else:
            current_exponent = exps[current_exp_ind][0]
            if (i > current_exponent):
                vec.append(0)
            else: # i == current_exponent
                vec.append(d[exps[current_exp_ind]])
                current_exp_ind += 1
    
    vec.reverse()
    return Matrix(vec)

def lin_rep(P, Q, p): # Returns a linear representation of ct[P^nQ] mod p
    R.<t> = LaurentPolynomialRing(GF(p), 1)
    max_deg = max(R(P).degree() - 1, R(R(P)(t^(-1))).degree() - 1, R(Q).degree(), R(R(Q)(t^(-1))).degree())
    
    mats = []
    P_power = R(1)
    for k in range(p):
        mats.append(compute_mat_for_poly(P_power, max_deg, p))
        P_power = P_power*R(P)
    
    return (poly_to_vec(Q, max_deg), mats, poly_to_vec(R(1), max_deg).transpose())

def apply_lin_rep(v, mats, w, n, p): # Returns v * phi((n)_p) * w
    vec = v
    digits = Integer(n).digits(p)
    for digit in digits:
        vec = vec * mats[digit]
    return (vec * w)[0][0]

def lin_rep_to_machine(v, mats, p, state_bound): # Converts a linear representation to a DFA
    k = 0
    states = [v]
    transitions = [[]]

    while (k < len(states)):
        for i in range(p):
            next_state = states[k]*mats[i]
            transitions[k].append(next_state)
            
            found = false
            for j in range(len(states)):
                if (states[j] == next_state):
                    found = true
                    break
            if not found:
                states.append(next_state)
                transitions.append([])
                if (len(states) > state_bound):
                    raise Exception("Number of states exceeded " + str(state_bound))
        k = k+1
    
    return (states, transitions)

def serialize_lin_rep_machine(machine, w, p): # Serializes a machine (as output by lin_rep_to_machine) over alphabet GF(p) for use with Walnut
    (states, transitions) = machine
    s = "lsd_" + str(p) + "\n\n"
    for k in range(len(states)):
        s += str(k) + " " + str((states[k]*w)[0][0]) + "\n"
        for i in range(len(transitions[k])):
            next_state = transitions[k][i]
            s += str(i) + " -> " + str(states.index(next_state)) + "\n"
        s += "\n"
    return s

def compute_shortest_prop(P, Q, p, prop, state_bound): # Returns the first value for which ct[P^nQ] satisfies prop
    # It is asymptotically faster to build the DFA than to loop through values of n
    # A new implementation is required since this library constructs lsd DFAs
    # but finding the first prop-satisfying value requires minimizing with a msd DFA.
    if (prop(Q.constant_coefficient())):
        return 0
    
    (v, mats, w) = lin_rep(P, Q, p)
    k = 0
    states = [(w, [])]
    transitions = [[]]

    while (k < len(states)):
        for i in range(p):
            (state, path) = states[k]
            next_state = mats[i]*state

            if (prop((v*next_state)[0][0])):
                new_path = path.copy()
                new_path.append(i)
                new_path.reverse()
                return Integer(new_path, p)

            transitions[k].append(next_state)
            
            found = false
            for j in range(len(states)):
                (previous_state, _) = states[j]
                if (previous_state == next_state):
                    found = true
                    break
            if not found:
                new_path = path.copy()
                new_path.append(i)
                states.append((next_state, new_path))
                transitions.append([])
                if (len(states) > state_bound):
                    raise Exception("Number of states exceeded " + str(state_bound))
        k = k+1
    
    return None

def compute_shortest_element(P, Q, p, element, state_bound): # Returns the first value for which p divides element - ct[P^nQ]
    def equals_element(n):
        return (n == element)
    
    return compute_shortest_prop(P, Q, p, equals_element, state_bound)

def compute_shortest_zero(P, Q, p, state_bound): # Returns the first value of n for which p divides ct[P^nQ]
    import DFA
    try: # First optimistically try to compute the lsd-first DFA with a low state_bound
        if (DFA.PolyAutoFailOnZero(P, Q, p, 100+p) != None):
            return None
    except:
        None # Do nothing here
    
    return compute_shortest_element(P, Q, p, 0, state_bound)

def compute_shortest_non_zero(P, Q, p, state_bound): # Returns the first value of n for which p does not divide ct[P^nQ]]
    def not_zero(n):
        return (n != 0)
    
    return compute_shortest_prop(P, Q, p, not_zero, state_bound)