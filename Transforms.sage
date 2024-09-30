from PolyUtil import compute_triangle

def transform(seq, P, p, triangle = None, offset = None): # Constructs a sequence with generating function Q such that seq[n] = ct[P^n*Q] mod p
    # P must be of the form t^-1 + a_0*1 + ... + a_r*t^r (with r >= 0)

    F = GF(p)
    transformed = [F(seq[0])]

    if (triangle == None):
        (triangle, offset) = compute_triangle(P, p, len(seq))
    
    for i in range(1, len(seq)):
        sum = F(seq[i])
        for j in range(i):
            sum -= triangle[i][offset - j]*transformed[j]
        transformed.append(F(sum))
    return transformed

def transform_inverse(seq, P, p, triangle = None, offset = None): # Treating seq as the coefficients of Q, returns the sequence ct[P^nQ] mod p
    F = GF(p)
    untransformed = [F(seq[0])]

    if (triangle == None):
        (triangle, offset) = compute_triangle(P, p, len(seq))
    
    for i in range(1, len(seq)):
        sum = F(0)
        for j in range(i+1):
            sum += triangle[i][offset - j]*seq[j]
        untransformed.append(F(sum))
    return untransformed

def seq_from_string(num_str): # Turns each character of num_str as an element of an integer sequence
    nums = []
    for c in num_str:
        nums.append(int(c))
    return nums

def DFA_guess(seq, P, p): # Given a long enough prefix, seq, returns a DFA for seq if seq is automatic
    # P must be of the form t^-1 + a_0*1 + ... + a_r*t^r (with r >= 0)
    # E.g. P = t^-1 + 1 will work!

    # The algorithm is to transform seq under P so that seq[i] = ct[P^i*Q] mod p
    # where Q is the P-transform of seq, then a_{pn + k} = ct[P^n*Lambda_p(P^kQ)] mod p so
    # that we can build a DFA from the relation Q -k-> Lambda_p(P^kQ) with output Q -> ct[Q]

    F = GF(p)
    (triangle, offset) = compute_triangle(P, p, max(len(seq), p))

    def possibly_equal(seq1, seq2): # True if one is a prefix of the other
        for i in range(min(len(seq1), len(seq2))):
            if (seq1[i] != seq2[i]):
                return false
        return true

    def step(ds, k): # Computes Lambda_p(P^k*ds)
        result = []
        for i in range(len(ds)//p - 1): # Compute the coefficient of t^(p*i) in P^k*ds
            sum = F(0)
            for j in range(-k, k*P.degree() + 1):
                if (p*i - j >= 0):
                    sum += triangle[k][offset + j]*ds[p*i - j]
            result.append(F(sum))
        return result
    
    states = [transform(seq, P, p, triangle, offset)] # The transforms of the p-kernel sequences
    transitions = []
    current_index = 0

    while (current_index < len(states)):
        current_state = states[current_index]
        new_transitions = []
        for i in range(p):
            new_state = step(current_state, i)

            new_state_index = len(states)
            for j in range(len(states)): # See if this state has already been found
                if (possibly_equal(states[j], new_state)):
                    new_state_index = j
                    break
            if (new_state_index == len(states)): # If it is new, add it to states
                states.append(new_state)
            
            new_transitions.append((i, new_state_index))

        transitions.append(new_transitions)
        current_index += 1
    
    def output_func(state):
        return state[0] # state[0] is the constant term
    
    return (states, transitions, output_func)