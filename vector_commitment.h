//
// Created by root on 23-4-25.
//

#ifndef NTL_TEST_VECTOR_COMMITMENT_H
#define NTL_TEST_VECTOR_COMMITMENT_H
#include "accumulator.h"

class vector_commitment {
public:
    struct VectorProof {
        accumulator::MembershipProof membership_proof;
        accumulator::NonmembershipProof nonmembership_proof;
    };
    pair<vector<mpz_class>, vector<mpz_class>> group_elems_by_bit(vector<pair<bool, mpz_class>> bits);
    void empty(accumulator::Accumulator &self);
    pair<accumulator::Accumulator, VectorProof> update(accumulator::Accumulator vc, vector<mpz_class> vc_acc_set, vector<pair<bool, mpz_class>> bits);
    VectorProof  open(accumulator::Accumulator vc, vector<mpz_class> vc_acc_set, vector<mpz_class> zero_bits, vector<pair<mpz_class , accumulator::Witness>> one_bit_witnesses);
    bool verify(accumulator::Accumulator vc, vector<pair<bool, mpz_class>> bits, VectorProof vc_proof);
};

#endif //NTL_TEST_VECTOR_COMMITMENT_H