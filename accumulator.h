//
// Created by root on 23-4-25.
//

#ifndef NTL_TEST_ACCUMULATOR_H
#define NTL_TEST_ACCUMULATOR_H
#include "group/classgroup.h"
#include "proof/poe.h"
#include "proof/poke2.h"

class accumulator {
public:
    struct Accumulator{
        classgroup::classElem value;
        mpz_class phantom;
    };
    struct Witness{
        Accumulator witness;
    };
    struct MembershipProof{
            Witness witness;
            poe::Poe proof;
    };
    struct NonmembershipProof {
        mpz_class phantom;
        classgroup::classElem d, v, gv_inv;
        poke2::Poke2 poke2_proof;
        poe::Poe poe_proof;
    };
    bool eq(accumulator::Accumulator acc1, accumulator::Accumulator acc2);
    void empty(accumulator::Accumulator &acc);
    accumulator::Accumulator add(Accumulator acc, vector<mpz_class> elems);
    pair<accumulator::Accumulator,accumulator::MembershipProof> add_with_proof(Accumulator acc, vector<mpz_class> elems);
    accumulator::Accumulator divide_and_conquer(vector<pair<mpz_class, accumulator::Witness>> elem_witnesses);
    accumulator::Accumulator Delete(Accumulator acc, vector<pair<mpz_class , Witness>> elem_witnesses);
    pair<accumulator::Accumulator,accumulator::MembershipProof> delete_with_proof(Accumulator acc, vector<pair<mpz_class , Witness>> elem_witnesses);
    accumulator::MembershipProof prove_membership(Accumulator acc, vector<pair<mpz_class , Witness>> elem_witnesses);
    bool verify_membership(Accumulator acc, mpz_class t, MembershipProof mbp);
    bool verify_membership_batch(Accumulator acc, vector<mpz_class> elems, MembershipProof mbp);
    accumulator::Witness update_membership_witness(Accumulator acc, Witness witness, vector<mpz_class> tracked_elems, vector<mpz_class> untracked_additions, vector<mpz_class> untracked_deletions);
    accumulator::NonmembershipProof prove_nonmembership(Accumulator acc, vector<mpz_class> acc_set, vector<mpz_class> elems);
    bool verify_nonmembership(Accumulator acc, vector<mpz_class> elems, NonmembershipProof nmbp);
    accumulator::Accumulator new_acc(vector<mpz_class> data);
    accumulator::Witness compute_subset_witness(Witness wit, vector<mpz_class> witness_set, vector<mpz_class> witness_subset);
    vector<pair<mpz_class, accumulator::Witness>> compute_individual_witnesses(Witness wit, vector<mpz_class> elems);
    vector<pair<mpz_class, accumulator::Witness>> root_factor(Witness wit, vector<mpz_class> elems);
    void test_add();
    void test_delete();
    void test_delete_empty();
    void test_delete_bad_witness();
    void test_update_membership_witness();
    void test_update_membership_witness_failure();
    void test_prove_nonmembership();
    void test_compute_sub_witness();
    void test_compute_sub_witness_failure();
    void test_compute_individual_witnesses();
private:
    pair<accumulator::Accumulator, mpz_class> add_(Accumulator acc, vector<mpz_class> elems);
    pair<accumulator::Accumulator, mpz_class> delete_(Accumulator acc, vector<pair<mpz_class , Witness>> elem_witnesses);
};


#endif //NTL_TEST_ACCUMULATOR_H
