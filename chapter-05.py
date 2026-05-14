import random
from functools import reduce

from py_ecc.bn128 import FQ, Z1, add, eq, multiply
from py_ecc.bn128 import curve_order as p


def rand():
    return random.randint(0, p)


def add_points(*points):
    return reduce(add, points, Z1)


# if points = G1, G2, G3, G4 and scalars = a,b,c,d vector_commit returns
# aG1 + bG2 + cG3 + dG4
def vector_commit(points, scalars):
    return reduce(add, [multiply(P, i) for P, i in zip(points, scalars)], Z1)


def inner(a, b):
    res = 0
    for ai, bi in zip(a, b):
        res += ai * bi
    return res % p


# these EC points have unknown discrete logs:
G = [
    (
        FQ(
            6286155310766333871795042970372566906087502116590250812133967451320632869759
        ),
        FQ(
            2167390362195738854837661032213065766665495464946848931705307210578191331138
        ),
    ),
    (
        FQ(
            6981010364086016896956769942642952706715308592529989685498391604818592148727
        ),
        FQ(
            8391728260743032188974275148610213338920590040698592463908691408719331517047
        ),
    ),
    (
        FQ(
            15884001095869889564203381122824453959747209506336645297496580404216889561240
        ),
        FQ(
            14397810633193722880623034635043699457129665948506123809325193598213289127838
        ),
    ),
    (
        FQ(
            6756792584920245352684519836070422133746350830019496743562729072905353421352
        ),
        FQ(
            3439606165356845334365677247963536173939840949797525638557303009070611741415
        ),
    ),
]

H = [
    (
        FQ(
            13728162449721098615672844430261112538072166300311022796820929618959450231493
        ),
        FQ(
            12153831869428634344429877091952509453770659237731690203490954547715195222919
        ),
    ),
    (
        FQ(
            17471368056527239558513938898018115153923978020864896155502359766132274520000
        ),
        FQ(
            4119036649831316606545646423655922855925839689145200049841234351186746829602
        ),
    ),
    (
        FQ(
            8730867317615040501447514540731627986093652356953339319572790273814347116534
        ),
        FQ(
            14893717982647482203420298569283769907955720318948910457352917488298566832491
        ),
    ),
    (
        FQ(419294495583131907906527833396935901898733653748716080944177732964425683442),
        FQ(
            14467906227467164575975695599962977164932514254303603096093942297417329342836
        ),
    ),
]

F = G[0]  # another point for v

B = (
    FQ(12848606535045587128788889317230751518392478691112375569775390095112330602489),
    FQ(18818936887558347291494629972517132071247847502517774285883500818572856935411),
)

# scalar multiplication example: multiply(G, 42)
# EC addition example: add(multiply(G, 42), multiply(G, 100))


# remember to do all arithmetic modulo p
def commit(a, sL, b, sR, alpha, beta, gamma, tau_1, tau_2):

    A = add(add_points(vector_commit(G, a)), add_points(vector_commit(H, b)))
    A = add(A, multiply(B, alpha))
    S = add(add_points(vector_commit(G, sL)), add_points(vector_commit(H, sR)))
    S = add(S, multiply(B, beta))
    V = add(multiply(F, (inner(a, b))), multiply(B, gamma))
    T1 = add(multiply(F, (inner(a, sR) + inner(b, sL))), multiply(B, tau_1))
    T2 = add(multiply(F, inner(sL, sR)), multiply(B, tau_2))
    return A, S, V, T1, T2


def evaluate(f_0, f_1, f_2, u):

    if type(f_0) is list:
        res = []
        for l in range(len(f_0)):
            res.append((f_0[l] + f_1[l] * u + f_2 * u**2) % p)

        return res

    elif type(f_0) is int:
        return (f_0 + f_1 * u + f_2 * u**2) % p
    else:
        return -1


def prove(blinding_0, blinding_1, blinding_2, u):
    pi = blinding_0 + blinding_1 * u + blinding_2 * u**2
    return pi % p


## step 0: Prover and verifier agree on G and B

## step 1: Prover creates the commitments
a = [14, 12, 21, 23]
b = [32, 14, 991, 24]
sL = [rand() for i in range(4)]
sR = [rand() for i in range(4)]
t1 = inner(a, sR) + inner(b, sL)
t2 = inner(sR, sL)
### blinding terms
alpha = rand()
beta = rand()
gamma = rand()
tau_1 = rand()
tau_2 = rand()

A, S, V, T1, T2 = commit(a, sL, b, sR, alpha, beta, gamma, tau_1, tau_2)

## step 2: Verifier picks u
u = rand()

## step 3: Prover evaluates l(u), r(u), t(u) and creates evaluation proofs
l_u = evaluate(a, sL, 0, u)
r_u = evaluate(b, sR, 0, u)
t_u = evaluate(inner(a, b), t1, t2, u)

pi_lr = prove(alpha, beta, 0, u)
pi_t = prove(gamma, tau_1, tau_2, u)

## step 4: Verifier accepts or rejects
assert t_u == inner(l_u, r_u), "tu !=<lu, ru>"
assert eq(
    add(A, multiply(S, u)),
    add_points(vector_commit(G, l_u), vector_commit(H, r_u), multiply(B, pi_lr)),
), "l_u or r_u not evaluated correctly"
assert eq(
    add(multiply(F, t_u), multiply(B, pi_t)),
    add_points(V, multiply(T1, u), multiply(T2, u**2 % p)),
), "t_u not evaluated correctly"
print("OK")
