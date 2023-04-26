//
// Created by root on 23-4-18.
//

#ifndef NTL_TEST_UTIL_H
#define NTL_TEST_UTIL_H
#include "group/classgroup.h"

class util {
public:
//    mpz_class prime_hash_product();
    std::optional<classgroup::classElem> shamir_trick(classgroup::classElem xth_root, classgroup::classElem yth_root, mpz_class x, mpz_class y);
    std::optional<std::pair<mpz_class, mpz_class>> solve_linear_congruence(mpz_class a, mpz_class b, mpz_class m);
    void test_linear_congruence_solver();
    void test_linear_congruence_solver_no_solution();
    void test_shamir_trick();
    void test_shamir_trick_failure();
};


#endif //NTL_TEST_UTIL_H
