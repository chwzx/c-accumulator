//
// Created by root on 23-4-21.
//
#include "poe.h"


poe::Poe poe::prove(classgroup::classElem base, mpz_class e, classgroup::classElem result) {
    //    mpz_class l = hash_to_prime(&(base, exp, result));
    mpz_class l = 1;
    mpz_class q = e / l;
    classgroup cl;
    return {cl.exp(base, q)};
}

bool poe::verify(classgroup::classElem base, mpz_class e, classgroup::classElem result, poe::Poe proof) {
    //    mpz_class l = hash_to_prime(&(base, exp, result));
    mpz_class l = 1;
    mpz_class r;
    // e = ql + r
    r = e % l;
    // w = Q^l * u^r (Q = g ^ q)
    classgroup cl;
    classgroup::classElem w = cl.op(cl.exp(proof.Q, l), cl.exp(base, r));
    return cl.eq(w, result);
}

void poe::test_poe_small_exp() {
    // result = (g^2)^20 = g^40
    classgroup cl;
    classgroup::classElem base = cl.unknown_order_elem();
    // base = g^2
    base = cl.exp(base, 2);
    // g^40
    mpz_class exp = 20;
    classgroup::classElem result = cl.exp(base, exp);
    poe::Poe proof = poe::prove(base, exp, result);
    assert(poe::verify(base, exp, result, proof));
    assert(cl.eq(proof.Q, cl.exp(base, exp)));

    // (g^2)^35 = g^70
    mpz_class exp_2 = 35;
    classgroup::classElem result_2 = cl.exp(base, 35);
    poe::Poe proof_2 = poe::prove(base, exp_2, result_2);
    assert(poe::verify(base, exp_2, result_2, proof_2));
    assert(cl.eq(proof_2.Q, cl.exp(base, exp_2)));
}
