from PolyUtil import Lambda

def PolyAutoGeneral(P, Q, p, state_bound, fail_on_zero): # Returns state machine for ct[P^nQ] mod p as a list of states labeled by polynomials Q' and correspondingly indexed transitions
    R.<t> = LaurentPolynomialRing(GF(p), 1)
    
    if(fail_on_zero and R(Q).constant_coefficient() == 0):
        return None
    
    k = 0 # The index of the current state being processed for transitions
    states = [R(Q)] # Begin with Q state
    transitions = [[]] # Begin with no transitions
    
    while (k < len(states)): # Continue while there are new states Q' = states[k]
        for i in range(p): # Add transitions for states[k]
            new_state = Lambda(R(P)^i*states[k], p)
            
            if (fail_on_zero and R(simp).constant_coefficient() == 0): # If a zero state is found, return None if fail_on_zero is true
                    return None
            
            new_state_index = len(states) # Determine if new_state is a state we've already seen
            for j in range(len(states)):
                if (states[j] == new_state):
                    new_state_index = j
                    break
            if (new_state_index == len(states)): # If new_state is a new state, add it to the state list and abort of state_bound is broken
                states.append(new_state)
                transitions.append([])
                if (len(states) > state_bound):
                    raise Exception("Number of states exceeded " + str(state_bound))
            
            transitions[k].append((i, new_state_index))
        k = k+1

    def output_func(state):
        return R(state).constant_coefficient()
    
    return (states, transitions, output_func)

def PolyAuto(P, Q, p, state_bound): # Returns state machine for ct[P^nQ] mod p as a list of states labeled by polynomials Q' and correspondingly indexed transitions
    return PolyAutoGeneral(P, Q, p, state_bound, False)

def PolyAutoFailOnZero(P, Q, p, state_bound): # Returns None if a zero state is found. Otherwise the same as def PolyAuto above
    return PolyAutoGeneral(P, Q, p, state_bound, True)

def serialize(machine, p): # Serializes a machine over alphabet GF(p) for use with Walnut
    (states, transitions, output_func) = machine
    s = "lsd_" + str(p) + "\n\n"
    for k in range(len(states)):
        s += str(k) + " " + str(output_func(states[k])) + "\n"
        for i in range(len(transitions[k])):
            (j, ind) = transitions[k][i]
            s += str(j) + " -> " + str(ind) + "\n"
        s += "\n"
    return s

def evaluate(machine, p, input): # Evaluates machine on input integer
    digits = Integer(input).digits(p)
    (states, transitions, ouptut_func) = machine
    current_index = 0
    for digit in digits:
        (i, j) = transitions[current_index][digit]
        assert(i == digit)
        current_index = j
    return output_func(states[current_index])

# TODO: Add multi-variate support (not urgent)
# TODO: Tests in another file with some framework or something
# TODO: Benchmarks for functions computing equal values
# TODO: Make nice tutorial notebook(s)