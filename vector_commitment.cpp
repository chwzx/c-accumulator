//
// Created by root on 23-4-25.
//

#include "vector_commitment.h"
#include <set>

void vector_commitment::empty(accumulator::Accumulator &self) {
    accumulator ac;
    ac.empty(self);
}

pair<vector<mpz_class>, vector<mpz_class>> vector_commitment::group_elems_by_bit(vector<pair<bool, mpz_class>> bits) {
    vector<mpz_class> elems_with_one;
    vector<mpz_class> elems_with_zero;
    set<mpz_class> seen_indices;
    for(auto x : bits){
        if(seen_indices.count(x.second))
//            assert(0);
            return {};
        if(x.first){
            elems_with_one.push_back(x.second);
        } else {
            elems_with_zero.push_back(x.second);
        }
    }
    return {elems_with_zero, elems_with_one};
}

pair<accumulator::Accumulator, vector_commitment::VectorProof>
vector_commitment::update(accumulator::Accumulator vc, vector<mpz_class> vc_acc_set,
                          vector<pair<bool, mpz_class>> bits) {
    pair<vector<mpz_class>, vector<mpz_class>> z_o = group_elems_by_bit(bits);
    vector<mpz_class> elems_with_zero = z_o.first;
    vector<mpz_class> elems_with_one = z_o.second;
    accumulator acc;
    pair<accumulator::Accumulator,accumulator::MembershipProof> v_add_res = acc.add_with_proof(vc, elems_with_one);
    accumulator::Accumulator new_acc = v_add_res.first;
    accumulator::MembershipProof membership_proof = v_add_res.second;
    accumulator::NonmembershipProof nonmembership_proof = acc.prove_nonmembership(new_acc, vc_acc_set, elems_with_zero);
    return {new_acc, {membership_proof, nonmembership_proof}};
}

vector_commitment::VectorProof
vector_commitment::open(accumulator::Accumulator vc, vector<mpz_class> vc_acc_set, vector<mpz_class> zero_bits,
                        vector<pair<mpz_class, accumulator::Witness>> one_bit_witnesses) {
    accumulator acc;
    accumulator::MembershipProof membership_proof = acc.prove_membership(vc, one_bit_witnesses);
    accumulator::NonmembershipProof nonmembership_proof = acc.prove_nonmembership(vc, vc_acc_set, zero_bits);

    return {membership_proof, nonmembership_proof};
}

bool vector_commitment::verify(accumulator::Accumulator vc, vector<pair<bool, mpz_class>> bits,
                               vector_commitment::VectorProof vc_proof) {
    pair<vector<mpz_class>, vector<mpz_class>> group_result = group_elems_by_bit(bits);
    vector<mpz_class> elems_with_zero = group_result.first;
    vector<mpz_class> elems_with_one = group_result.second;
    accumulator acc;
    return acc.verify_membership_batch(vc, elems_with_one, vc_proof.membership_proof) && acc.verify_nonmembership(vc, elems_with_zero, vc_proof.nonmembership_proof);
}
