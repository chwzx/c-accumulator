//
// Created by root on 23-4-24.
//

#include "poke2.h"

poke2::Poke2 poke2::prove(classgroup::classElem base, mpz_class e, classgroup::classElem result) {
    classgroup cl;
    classgroup::classElem g = cl.unknown_order_elem();
    classgroup::classElem z = cl.exp(g, e);
//    let l = hash_to_prime(&(base, result, &z));
    mpz_class l = 1;
//    let alpha = blake2b(&(base, result, &z, &l));
    mpz_class alpha = 1;
    mpz_class q, r;
    mpz_divmod(q.get_mpz_t(), r.get_mpz_t(), e.get_mpz_t(), l.get_mpz_t());
    classgroup::classElem Q = cl.exp(cl.op(base, cl.exp(g, alpha)), q);
    return { z, Q, r };
}

bool poke2::verify(classgroup::classElem base, classgroup::classElem result, poke2::Poke2 proof) {
    classgroup cl;
    classgroup::classElem g = cl.unknown_order_elem();
//    let l = hash_to_prime(&(base, result, &z));
//    let alpha = blake2b(&(base, result, &z, &l));
    mpz_class l = 1;
    mpz_class alpha = 1;
    classgroup::classElem lhs = cl.op(
            cl.exp(proof.Q, l),
            cl.exp(cl.op(base, cl.exp(g, alpha)), proof.r)
    );
    classgroup::classElem rhs = cl.op(result, cl.exp(proof.z, alpha));
    return cl.eq(lhs, rhs);
}

bool poke2::eq(poke2::Poke2 p1, poke2::Poke2 p2) {
    classgroup cl;
    return p1.r == p2.r && cl.eq(p1.z, p2.z) && cl.eq(p1.Q, p2.Q);
}



void poke2::test_poke2() {
    // 2^20 = 1048576
    classgroup cl;
    classgroup::classElem base = cl.unknown_order_elem();
    base = cl.exp(base, 2);
    mpz_class exp = 20;
    classgroup::classElem result = cl.exp(base, exp);
    poke2::Poke2 proof = poke2::prove(base, exp, result);
    assert(poke2::verify(base, result, proof));
    assert(eq(proof,{
                    z: cl.exp(base, 10),
                    Q: cl.exp(base, 30),
                    r: 0
        }
    ));

    // 2^35 = 34359738368
    mpz_class exp_2 = 35;
    classgroup::classElem result_2 = cl.exp(base, 35);
    poke2::Poke2 proof_2 = poke2::prove(base, exp_2, result_2);
    assert(poke2::verify(base, result_2, proof_2));
    // Cannot verify wrong base/exp/result triple with wrong pair.
    assert(!poke2::verify(base, result_2, proof));
    assert(eq(proof_2,{
                      z: cl.exp(cl.unknown_order_elem(), 35),
                      Q: cl.exp(cl.unknown_order_elem(), 105),
                      r: 0
              }
    ));
}


void poke2::test_poke2_negative() {

}



