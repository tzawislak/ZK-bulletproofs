import random

from py_ecc.bn128.bn128_curve import FQ, add, multiply
from py_ecc.bn128.bn128_curve import curve_order as p


def random_field_element():
    return random.randint(0, p)


# these EC points have unknown discrete logs:
G = (
    FQ(6286155310766333871795042970372566906087502116590250812133967451320632869759),
    FQ(2167390362195738854837661032213065766665495464946848931705307210578191331138),
)


B = (
    FQ(12848606535045587128788889317230751518392478691112375569775390095112330602489),
    FQ(18818936887558347291494629972517132071247847502517774285883500818572856935411),
)

# scalar multiplication example: multiply(G, 42)
# EC addition example: add(multiply(G, 42), multiply(G, 100))

# remember to do all arithmetic modulo p


def commit(f_0, f_1, f_2, gamma_0, gamma_1, gamma_2, G, B):
    print(f_0)
    C0 = add(multiply(G, f_0), multiply(B, gamma_0))
    C1 = add(multiply(G, f_1), multiply(B, gamma_1))
    C2 = add(multiply(G, f_2), multiply(B, gamma_2))

    # return the commitments as a tuple (C0, C1, C2)
    return C0, C1, C2


def evaluate(f_0, f_1, f_2, u):
    return (f_0 + f_1 * u + f_2 * u**2) % p


def prove(gamma_0, gamma_1, gamma_2, u):
    # fill this in
    pi = (gamma_0 + gamma_1 * u + gamma_2 * u**2) % p
    return pi


def verify(C0, C1, C2, G, B, f_u, pi, y):
    # Return true or false
    c0puc1 = add(C0, multiply(C1, f_u))
    committed = add(c0puc1, multiply(C2, f_u**2))
    challenge = add(multiply(G, y), multiply(B, pi))

    return committed == challenge


# step 0: Prover and verifier agree on G and B

# step 1: Prover creates the commitments
# f(x) = f_0 + f_1x + f_2x^2
f_0 = 1
f_1 = 2
f_2 = 3

# blinding terms
gamma_0 = random_field_element()
gamma_1 = random_field_element()
gamma_2 = random_field_element()

C0, C1, C2 = commit(f_0, f_1, f_2, gamma_0, gamma_1, gamma_2, G, B)

# step 2: Verifier picks u
u = random_field_element()

# step 3: Prover evaluates f(u) and pi

y = evaluate(f_0, f_1, f_2, u)
pi = prove(gamma_0, gamma_1, gamma_2, u)

# step 4: Verifier accepts or rejects
if verify(C0, C1, C2, G, B, u, pi, y):
    print("accept")
else:
    print("reject")
