def main():
    if len(sys.argv) < 4:
        print("Usage: sage inverse_polynomial.sage <coeffs_a> <coeffs_phi> <q>")
        return
    f = int(sys.argv[2])
    phi = euler_phi(f)
    q = int(sys.argv[3])
    coeffs_a = [int(c) for c in sys.argv[1].strip('[]').split(',')]
    R = PolynomialRing(GF(q), 'x')
    P = R.cyclotomic_polynomial(f);
    poly_a = R(coeffs_a)
    pari_poly_a = pari(poly_a)
    pari_poly_p = pari(P)
    inv_c_q_pari = pari(1 / pari_poly_a).Mod(pari_poly_p).lift()

    if inv_c_q_pari is None:
        print("None")
    else:
        inv = []
        for v in inv_c_q_pari:
            inv.append(int(v))
        print(inv)

if __name__ == "__main__":
    main()

