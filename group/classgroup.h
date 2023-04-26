//
// Created by root on 23-4-18.
//

#ifndef NTL_TEST_CLASSGROUP_H
#define NTL_TEST_CLASSGROUP_H
#include <gmpxx.h>
#include <iostream>
#include <cassert>
#include <gmp.h>

using namespace std;

static const char* CLASS_GROUP_DISCRIMINANT = "-30616069034807523947093657516320815215492876376165067902716988657802400037331914448218251590830\
  1102189519215849430413184776658192481976276720778009261808832630304841711366872161223643645001916\
  6969493423497224870506311710491233557329479816457723381368788734079933165653042145718668727765268\
  0575673207678516369650123480826989387975548598309959486361425021860161020248607833276306314923730\
  9854570972702350567411779734372573754840570138310317754359137013512655926325773048926718050691092\
  9453371727344087286361426404588335160385998280988603297435639020911295652025967761702701701471162\
  3966286152805654229445219531956098223";

static const mpz_class DETAL = mpz_class (CLASS_GROUP_DISCRIMINANT);

class classgroup {
public:
    struct classElem{
        mpz_class a, b, c;
    };
    classElem normalize(mpz_class a, mpz_class b, mpz_class c);
    classElem reduce(mpz_class a, mpz_class b, mpz_class c);
    classElem square(classElem x); // no
    mpz_class discriminant(mpz_class a, mpz_class b, mpz_class c);
    bool validate(mpz_class a, mpz_class b, mpz_class c);
    bool is_reduced(mpz_class a, mpz_class b, mpz_class c);
    bool is_normal(mpz_class a, mpz_class b, mpz_class c);
    classElem op(classElem a, classElem b);
    classElem id();
    classElem inv(classElem a);
    classElem exp(classElem a, mpz_class n);
    classElem unknown_order_elem();
    bool eq(classElem self, classElem other);
    classElem elem(classElem abc);
    classElem construct_raw_elem_from_strings(char* a, char* b, char* c);
    void test_bad_elem();
    void test_elem_from();
    void test_equality();
    void test_reduce_basic();
    void test_normalize_basic();
    void test_discriminant_basic();
    void test_discriminant_across_ops();
    void test_op_single();
    void test_op_alternating();
    void test_op_complex();
    void test_id_basic();
    void test_id_repeated();
    void test_inv();
    void test_exp_basic();
    void test_square_basic();
private:
    classElem op_(classElem x, classElem y);
    classElem id_(mpz_class d);
    classElem inv_(classElem x);
    classElem exp_(classElem a, mpz_class n);
    classElem unknown_order_elem_(mpz_class d);
};


#endif //NTL_TEST_CLASSGROUP_H
