//
// Created by root on 23-4-21.
//

#ifndef NTL_TEST_POE_H
#define NTL_TEST_POE_H
#include "../group/classgroup.h"

class poe {
public:
    struct Poe{
        classgroup::classElem Q;
    };
    Poe prove(classgroup::classElem base, mpz_class e, classgroup::classElem result);
    bool verify(classgroup::classElem base, mpz_class e, classgroup::classElem result, Poe proof);
    void test_poe_small_exp();
};


#endif //NTL_TEST_POE_H
