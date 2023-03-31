import sys


def to_limbs(x, n):
    mask = 2 ^ 64 - 1
    return [(x >> (64 * i)) & mask for i in range(n)]


def curve_params():
    args = sys.argv[1:]
    if len(args) == 4:
        p, a, b = [int(x, 16) if x.startswith("0x") else int(x) for x in args[:-1]]
        filename = args[-1]
    elif len(args) > 0:
        print(
            f"""Usage for curve Y^2 = X^3 + aX + b over basefield GF(p). Integers in base 10 or 16 with 0x prefix.
    sage {sys.argv[0]} p a b output_filename

By default
    sage {sys.argv[0]}
generates parameters for the BN254 curve.
            """
        )
        sys.exit()
    else:
        # BN254 parameters
        p, a, b, filename = (
            0x30644E72E131A029B85045B68181585D97816A916871CA8D3C208C16D87CFD47,
            1,
            5612291247948481584627780310922020304781354847659642188369727566000581075360,
            "bn254",
        )
    return p, a, b, filename


p, a, b, filename = curve_params()
F = GF(p)
E = EllipticCurve(F, [a, b])
O = E.order()
n = O.p_primary_part(2)
log_n = n.log(2)
print(f"Curve's scalar field 2-adicity: {log_n}")

# Note that the algorithm below assume that E has a 2^log_n cyclic subgroup
# not just a 2^log_n subgroup
while True:
    g = E.random_point()
    G = (O // n) * g
    for i in range(log_n):
        if ((2 ^ i) * G).is_zero():
            print(f"  {2^i} * G is zero, trying a new point")
            break
    else:
        break

assert (n * G).is_zero()
print(f"found element G that generates the subgroup of order {2^log_n}")
R = E.random_element()
H = [R + i * G for i in range(2 ^ log_n)]
L = [h.xy()[0] for h in H]
S = [L[i] for i in range(0, n, 2)]
S_prime = [L[i] for i in range(1, n, 2)]

# Number of 64-bit limbs needed to represent field elements
num_limbs = (p.bit_length() + 63) // 64
print(f"number of limbs: {num_limbs}")

coset_file_name = f"{filename}_coset"
s = "\n".join([str(l) for x in L for l in to_limbs(int(x), num_limbs)])
with open(coset_file_name, "w") as f:
    f.write(s)

print(f"file {coset_file_name} written")


def isogenies(log_n, S, S_prime, E):
    isos = []
    for i in range(log_n, 0, -1):
        n = 1 << i
        nn = n // 2

        for iso in E.isogenies_prime_degree(2):
            psi = iso.x_rational_map()
            if len(set([psi(x) for x in S])) == nn:
                break
        isos.append(psi)
        S = [psi(x) for x in S[:nn]]
        S_prime = [psi(x) for x in S_prime[:nn]]
        E = iso.codomain()

    return isos


isos = isogenies(log_n - 1, S, S_prime, E)
s = "\n".join(
    [
        str(l)
        for psi in isos
        for coeff in list(psi.numerator()) + list(psi.denominator())
        for l in to_limbs(int(coeff), num_limbs)
    ]
)

isogenies_file_name = f"{filename}_isogenies"
with open(isogenies_file_name, "w") as f:
    f.write(s)
print(f"file {isogenies_file_name} written")
