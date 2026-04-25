from hashlib import sha256

from libnum import has_sqrtmod_prime_power, sqrtmod_prime_power
from py_ecc.bn128 import FQ, is_on_curve
from py_ecc.fields import field_properties

field_mod = field_properties["bn128"]["field_modulus"]


def generate_ec_point_basis(n: int):
    b = 3  # for bn128, y^2 = x^3 + 3
    seed = "RareSkills_seed"

    # we intend to generate n EC points that will serve the Pedersen commitment scheme
    # Step 1:
    # calculate a random field element

    if not hasattr(generate_ec_point_basis, "x"):
        generate_ec_point_basis.x = (
            int(sha256(seed.encode("ascii")).hexdigest(), 16) % field_mod
        )

    # rename as x for convenience
    x = generate_ec_point_basis.x
    entropy = 0

    vector_basis = []

    for nn in range(n):
        # Step 2: check if the generated field element belongs to the curve
        while not has_sqrtmod_prime_power((x**3 + b) % field_mod, field_mod, 1):
            # increment x, so hopefully we are on the curve
            x = (x + 1) % field_mod
            entropy = entropy + 1

        # pick the upper or lower point depending on if entropy is even or odd
        y = list(sqrtmod_prime_power((x**3 + b) % field_mod, field_mod, 1))[
            entropy & 1 == 0
        ]
        point = (FQ(x), FQ(y))
        assert is_on_curve(point, b), "sanity check"
        vector_basis.append(point)

    # new x value, for another basis generation
    generate_ec_point_basis.x = (
        int(sha256(str(x).encode("ascii")).hexdigest(), 16) % field_mod
    )
    print(f"Pedersen commitment basis: \n {vector_basis}")
    return vector_basis


vb = generate_ec_point_basis(1)
