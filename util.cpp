//
// Created by root on 23-4-18.
//

#include "util.h"
#include <iostream>
#include <cassert>

using namespace std;

std::optional<classgroup::classElem>
util::shamir_trick(classgroup::classElem xth_root, classgroup::classElem yth_root, mpz_class x, mpz_class y) {
    classgroup cl;
    if(!cl.eq(cl.exp(xth_root, x),cl.exp(yth_root, y))) {
        return {};
    }

    mpz_class gcd, a, b;
    mpz_gcdext(gcd.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t());
    if(gcd != 1){
        return {};
    }
    return {cl.op(cl.exp(xth_root, b), cl.exp(yth_root, a))};
}

std::optional<std::pair<mpz_class, mpz_class>> util::solve_linear_congruence(mpz_class a, mpz_class b, mpz_class m) {
    std::pair<mpz_class, mpz_class> res;
    mpz_class g, d, e;
    // g = gcd(a, m) => da + em = g
    mpz_gcdext(g.get_mpz_t(), d.get_mpz_t(), e.get_mpz_t(), a.get_mpz_t(), m.get_mpz_t());
    // q = floor_div(b, g)
    // r = b % g
    mpz_class q, r;
    mpz_divmod(q.get_mpz_t(), r.get_mpz_t(), b.get_mpz_t(), g.get_mpz_t());
    if(r != 0){
        return {};
    }

    mpz_class mu, v, tmp;
    // mu = qd % m
    mu = q * d % m;
    // v = m / g
    v = m / g;

    return {{mu,v}};
}

void util::test_linear_congruence_solver() {
    assert(solve_linear_congruence(3, 2, 4).value().first == -2 && solve_linear_congruence(3, 2, 4).value().second == 4);
    assert(solve_linear_congruence(5, 1, 2).value().first == 1 && solve_linear_congruence(5, 1, 2).value().second == 2);
    assert(solve_linear_congruence(2, 4, 5).value().first == -3 && solve_linear_congruence(2, 4, 5).value().second == 5);
    assert(solve_linear_congruence(230, 1081, 12167).value().first == 2491 && solve_linear_congruence(230, 1081, 12167).value().second == 529);
}

void util::test_linear_congruence_solver_no_solution() {
    assert(solve_linear_congruence(33, 7, 143) == nullopt);
    assert(solve_linear_congruence(13, 14, 39) == nullopt);
}

void util::test_shamir_trick() {
    mpz_class x = 13, y = 17, z = 19;
    classgroup cl;
    classgroup::classElem xth_root = cl.exp(cl.unknown_order_elem(), y * z);
    classgroup::classElem yth_root = cl.exp(cl.unknown_order_elem(), x * z);
    classgroup::classElem xyth_root = cl.exp(cl.unknown_order_elem(), z);
    if(util::shamir_trick(xth_root, yth_root, x, y) == nullopt)
        cout << "error\n";
    assert(cl.eq(util::shamir_trick(xth_root, yth_root, x, y).value(),xyth_root));
}

void util::test_shamir_trick_failure() {
    mpz_class x = 7, y = 14, z = 19;
    classgroup cl;
    classgroup::classElem xth_root = cl.exp(cl.unknown_order_elem(), y * z);
    classgroup::classElem yth_root = cl.exp(cl.unknown_order_elem(), x * z);
    assert(shamir_trick(xth_root, yth_root, x, y) == nullopt);
}

