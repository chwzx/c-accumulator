//
// Created by root on 23-4-24.
//

#ifndef NTL_TEST_POKE2_H
#define NTL_TEST_POKE2_H


#include "../group/classgroup.h"

class poke2 {
public:
    struct Poke2{
        classgroup::classElem z, Q;
        mpz_class r;
    };
    Poke2 prove(classgroup::classElem base, mpz_class e, classgroup::classElem result);
    bool verify(classgroup::classElem base, classgroup::classElem result, Poke2 proof);
    bool eq(Poke2 p1, Poke2 p2);
    void test_poke2();
    void test_poke2_negative();
};


#endif //NTL_TEST_POKE2_H
