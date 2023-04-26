//
// Created by root on 23-4-25.
//

#include "accumulator.h"
#include "util.h"

bool accumulator::eq(accumulator::Accumulator acc1, accumulator::Accumulator acc2) {
    classgroup cl;
    return cl.eq(acc1.value, acc2.value) && (acc1.phantom == acc2.phantom);
}

void accumulator::empty(accumulator::Accumulator &acc) {
    classgroup cl;
    acc.value = cl.unknown_order_elem();
    acc.phantom = 1;
}

pair<accumulator::Accumulator, mpz_class> accumulator::add_(accumulator::Accumulator acc, vector<mpz_class> elems) {
    classgroup cl;
//    let x = prime_hash_product(elems);
    mpz_class x = 1;
    for(auto n : elems)
        x *= n;
    acc.phantom *= x;
    acc.value = cl.exp(acc.value, x);
    // (g^old)^new  new
    return {acc, x};
}

accumulator::Accumulator accumulator::add(accumulator::Accumulator acc, vector<mpz_class> elems) {
    return add_(acc, elems).first;
}

pair<accumulator::Accumulator, accumulator::MembershipProof>
accumulator::add_with_proof(accumulator::Accumulator acc, vector<mpz_class> elems) {
    pair<accumulator::Accumulator, mpz_class> res = accumulator::add_(acc, elems);
    poe p;
    poe::Poe proof = p.prove(acc.value, res.second, res.first.value);
    return {res.first, {acc, proof}};
}

pair<accumulator::Accumulator, mpz_class>
accumulator::delete_(accumulator::Accumulator acc, vector<pair<mpz_class, accumulator::Witness>> elem_witnesses) {
    if(elem_witnesses.size() == 0)
        return {acc, 1};
//    let prime_witnesses = elem_witnesses.iter().map(|(elem, witness)| (hash_to_prime(elem), witness.0.value.clone())).collect::<Vec<_>>();
    classgroup cl;
    vector<pair<mpz_class, accumulator::Witness>> prime_witnesses = elem_witnesses;
    mpz_class prime_product;
    for(int i = 0; i < prime_witnesses.size(); i++){
        // membershipwitness
        if(!cl.eq(cl.exp(prime_witnesses[i].second.witness.value, prime_witnesses[i].first), acc.value)){
            assert(0);
            return {};
        }
    }
    Accumulator delete_elem_acc = divide_and_conquer(elem_witnesses);
    util ul;
    return {{ul.shamir_trick(delete_elem_acc.value, acc.value, delete_elem_acc.phantom, 1).value(), acc.phantom / delete_elem_acc.phantom},delete_elem_acc.phantom};
}

accumulator::Accumulator accumulator::divide_and_conquer(vector<pair<mpz_class, accumulator::Witness>> elem_witnesses) {
    if(elem_witnesses.size() == 1){
        return {elem_witnesses.front().second.witness.value, elem_witnesses.front().first};
    }
    int half = elem_witnesses.size() / 2;
    vector<pair<mpz_class, accumulator::Witness>> elem_witnesses_l, elem_witnesses_r;
    for(int i = 0; i < half; i++){
        elem_witnesses_l.push_back({elem_witnesses[i].first, elem_witnesses[i].second});
    }
    for(int i = half; i < elem_witnesses.size(); i++){
        elem_witnesses_l.push_back({elem_witnesses[i].first, elem_witnesses[i].second});
    }
    accumulator::Accumulator pv1 = divide_and_conquer(elem_witnesses_l);
    accumulator::Accumulator pv2 = divide_and_conquer(elem_witnesses_r);
    util u;
    return {u.shamir_trick(pv1.value, pv2.value, pv1.phantom, pv2.phantom).value(), pv1.phantom * pv2.phantom};
}

accumulator::Accumulator accumulator::Delete(accumulator::Accumulator acc, vector<pair<mpz_class, accumulator::Witness>> elem_witnesses) {
    return accumulator::delete_(acc, elem_witnesses).first;
}

pair<accumulator::Accumulator, accumulator::MembershipProof>
accumulator::delete_with_proof(accumulator::Accumulator acc,
                               vector<pair<mpz_class, accumulator::Witness>> elem_witnesses) {
    pair<accumulator::Accumulator, mpz_class> res = delete_(acc, elem_witnesses);
    poe p;
    poe::Poe proof = p.prove(res.first.value, res.second, acc.value);
    return {res.first, {res.first, proof}};
}

bool accumulator::verify_membership(accumulator::Accumulator acc, mpz_class t, accumulator::MembershipProof mbp) {
//    let exp = hash_to_prime(t);'
    mpz_class exp = t;
    poe p;
    return p.verify(mbp.witness.witness.value, exp, acc.value, mbp.proof);
}

accumulator::MembershipProof accumulator::prove_membership(accumulator::Accumulator acc,
                                                           vector<pair<mpz_class, accumulator::Witness>> elem_witnesses) {

    accumulator::Accumulator witness_accum = Delete(acc, elem_witnesses);
//    let prod = elem_witnesses.iter().map(|(t, _)| hash_to_prime(t)).product();
    mpz_class prod = 1;
    for(auto t : elem_witnesses){
        prod *= t.first;
    }
    poe p;
    poe::Poe proof = p.prove(witness_accum.value, prod, acc.value);
    return {witness_accum, proof};
}

bool accumulator::verify_membership_batch(accumulator::Accumulator acc, vector<mpz_class> elems,
                                          accumulator::MembershipProof mbp) {
//    let exp = prime_hash_product(elems);
    mpz_class exp = 1;
    for(auto elem : elems)
        exp *= elem;
    poe p;
    return p.verify(mbp.witness.witness.value, exp, acc.value, mbp.proof);
}

accumulator::Witness accumulator::update_membership_witness(accumulator::Accumulator acc, accumulator::Witness witness,
                                                            vector<mpz_class> tracked_elems,
                                                            vector<mpz_class> untracked_additions,
                                                            vector<mpz_class> untracked_deletions) {
//    mpz_class x = prime_hash_product(tracked_elems);
//    mpz_class x_hat = prime_hash_product(untracked_deletions);
    for( auto n : tracked_elems) {
        if(std::find(untracked_additions.begin(), untracked_additions.end(), n) != untracked_additions.end()
        || std::find(untracked_deletions.begin(), untracked_deletions.end(), n) != untracked_deletions.end()
        ) {
            assert(0);
            return {};
        }
    }
    mpz_class x = 1;
    for(auto n : tracked_elems)
        x *= n;

    mpz_class x_hat = 1;
    for(auto n : untracked_deletions)
        x_hat *= n;

    mpz_class gcd, a, b;
    mpz_gcdext(gcd.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t(), x.get_mpz_t(), x_hat.get_mpz_t());
    assert(gcd == 1);

    classgroup cl;
    accumulator::Witness w = {add(witness.witness, untracked_additions)};
    classgroup::classElem w_to_b = cl.exp(w.witness.value, b);
    classgroup::classElem acc_new_to_a = cl.exp(acc.value, a);
    return {{cl.op(w_to_b, acc_new_to_a), w.witness.phantom / x_hat}};
}

accumulator::NonmembershipProof
accumulator::prove_nonmembership(accumulator::Accumulator acc, vector<mpz_class> acc_set, vector<mpz_class> elems) {
//    mpz_class x = elems.iter().map(hash_to_prime).product();
    mpz_class x = 1;
    for(auto n : elems)
        x *= n;
//    mpz_class s = acc_set.iter().map(hash_to_prime).product();
    mpz_class s = 1;
    for(auto n : acc_set)
        s *= n;

    mpz_class gcd, a, b;
    mpz_gcdext(gcd.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t(), x.get_mpz_t(), s.get_mpz_t());
    assert(gcd == 1);

    classgroup cl;
    classgroup::classElem g = cl.unknown_order_elem();
    classgroup::classElem d = cl.exp(g, a);
    classgroup::classElem v = cl.exp(acc.value, b);
    classgroup::classElem gv_inv = cl.op(g, cl.inv(v));
    poke2 p2;
    poe p;
    poke2::Poke2 poke2_proof = p2.prove(acc.value, b, v);
    poe::Poe poe_proof = p.prove(d, x, gv_inv);
    return {1, d, v, gv_inv, poke2_proof, poe_proof};
}

bool accumulator::verify_nonmembership(accumulator::Accumulator acc, vector<mpz_class> elems,
                                       accumulator::NonmembershipProof nmbp) {
//    let x = elems.iter().map(hash_to_prime).product();
    mpz_class x = 1;
    for(auto n : elems)
        x *= n;
    poke2 p2;
    poe p;
    return p2.verify(acc.value, nmbp.v, nmbp.poke2_proof) && p.verify(nmbp.d, x, nmbp.gv_inv, nmbp.poe_proof);
}

accumulator::Accumulator accumulator::new_acc(vector<mpz_class> data) {
    accumulator::Accumulator acc;
    accumulator::empty(acc);
    return accumulator::add(acc, data);
}

void accumulator::test_add() {
    accumulator::Accumulator acc = new_acc({3, 5});
    vector<mpz_class> new_elems = {7, 11};
    pair<accumulator::Accumulator, accumulator::MembershipProof> res = add_with_proof(acc,new_elems);
    classgroup cl;
    classgroup::classElem acc_expected = cl.exp(
            cl.unknown_order_elem(),
            1155
    );
    assert(cl.eq(res.first.value, acc_expected));
    assert(accumulator::verify_membership_batch(res.first, new_elems, res.second));
}

void accumulator::test_delete() {
    Accumulator acc_0 = new_acc({3, 5});
    pair<accumulator::Accumulator, accumulator::MembershipProof> add_res = add_with_proof(acc_0, {7});
//    let (acc_1, c_proof) = acc_0.clone().add_with_proof(&["c"]);
    pair<accumulator::Accumulator, accumulator::MembershipProof> del_res = delete_with_proof(add_res.first, {{7, add_res.second.witness}});
//    let (acc_2, proof) = acc_1.clone().delete_with_proof(&[("c", c_proof.witness)]).expect("valid delete expected");
    assert(accumulator::eq(del_res.first, acc_0));
    assert(accumulator::verify_membership(add_res.first, 7, del_res.second));
}

void accumulator::test_delete_empty() {
    accumulator::Accumulator acc = new_acc({3, 5});
    pair<accumulator::Accumulator, accumulator::MembershipProof> del_res = delete_with_proof(acc, {});
    assert(accumulator::eq(del_res.first, acc));
    assert(accumulator::verify_membership_batch(acc, {}, del_res.second));
}

void accumulator::test_delete_bad_witness() {
    accumulator::Accumulator acc;
    empty(acc);
    Witness a_witness = {new_acc({5, 7})};
    Witness b_witness = {new_acc({3, 7})};
    Delete(acc, {{3, a_witness}, {5, b_witness}});
}

void accumulator::test_update_membership_witness() {
    accumulator::Accumulator acc = new_acc({3, 5, 7});
    cout << "acc:" << acc.phantom << endl;
    accumulator::Witness witness = {new_acc({7, 11})};
    accumulator::Witness witness_new = update_membership_witness(acc, witness, {3}, {5}, {11});
    cout << "witness_new:" << witness_new.witness.phantom << endl;
    assert(accumulator::eq(add(witness_new.witness, {3}),acc));
}

void accumulator::test_update_membership_witness_failure() {
    accumulator::Accumulator acc = new_acc({3, 5, 7});
    accumulator::Witness witness = {new_acc({7, 11})};
    update_membership_witness(acc, witness, {3}, {5}, {3});
}

void accumulator::test_prove_nonmembership() {
    vector<mpz_class> acc_set = {3, 5};
    accumulator::Accumulator acc = new_acc(acc_set);
    vector<mpz_class> non_members = {7, 11};

    accumulator::NonmembershipProof proof = prove_nonmembership(acc, acc_set, non_members);
    assert(accumulator::verify_nonmembership(acc, non_members, proof));
}

accumulator::Witness accumulator::compute_subset_witness(accumulator::Witness wit, vector<mpz_class> witness_set,
                                                         vector<mpz_class> witness_subset) {
    for(auto witness : witness_subset){
        if(std::find(witness_set.begin(), witness_set.end(), witness) == witness_set.end()){
            assert(0);
            return {};
        }
    }
//    let numerator = prime_hash_product(witness_set);
//    let denominator = prime_hash_product(witness_subset);
//    let (quotient, remainder) = numerator.div_rem(denominator);
    mpz_class numerator = 1;
    for(auto x : witness_set)
        numerator *= x;
    mpz_class denominator = 1;
    for(auto x : witness_subset)
        denominator *= x;

    mpz_class quotient, remainder;
    mpz_divmod(quotient.get_mpz_t(), remainder.get_mpz_t(), numerator.get_mpz_t(), denominator.get_mpz_t());
    assert(remainder == 0);
    classgroup cl;
    return {{cl.exp(wit.witness.value, quotient), numerator / denominator}};
}

vector<pair<mpz_class, accumulator::Witness>>
accumulator::compute_individual_witnesses(accumulator::Witness wit, vector<mpz_class> elems) {
//    mpz_class hashes = elems.iter().map(hash_to_prime).collect::<Vec<_>>();
    return root_factor(wit, elems);
}

vector<pair<mpz_class, accumulator::Witness>> accumulator::root_factor(accumulator::Witness wit, vector<mpz_class> elems) {
    if (elems.size() == 1) {
        return {{elems.front(), wit}};
    }
    int half_n = elems.size() / 2;
    classgroup cl;
    accumulator::Witness g_l;
    g_l = wit;
    vector<mpz_class> elems_l;
    for(int i = 0; i < half_n; i++){
        elems_l.push_back(elems[i]);
        g_l.witness.value = cl.exp(g_l.witness.value, elems[i]);
        g_l.witness.phantom *= elems[i];
    }

    accumulator::Witness g_r;
    g_r = wit;
    vector<mpz_class> elems_r;
    for(int i = half_n; i < elems.size(); i++){
        elems_r.push_back(elems[i]);
        g_r.witness.value = cl.exp(g_r.witness.value, elems[i]);
        g_r.witness.phantom *= elems[i];
    }

    vector<pair<mpz_class, accumulator::Witness>> L = root_factor(g_r, elems_l);
    vector<pair<mpz_class, accumulator::Witness>> R = root_factor(g_l, elems_r);
    L.insert(L.end(), R.begin(), R.end());
    return L;
}

void accumulator::test_compute_sub_witness() {
    accumulator::Witness empty_witness;
    empty(empty_witness.witness);
    accumulator::Witness sub_witness = compute_subset_witness(empty_witness, {3, 5}, {3});
    accumulator::Witness exp_quotient_expected = {new_acc({5})};
    assert(accumulator::eq(sub_witness.witness, exp_quotient_expected.witness));
}

void accumulator::test_compute_sub_witness_failure() {
    accumulator::Witness empty_witness;
    empty(empty_witness.witness);
    compute_subset_witness(empty_witness, {3, 5}, {11});
}

void accumulator::test_compute_individual_witnesses() {
    accumulator::Accumulator acc = new_acc({3, 5, 7, 11, 13});
    Witness witness_multiple = Witness({new_acc({3, 11})});
    vector<pair<mpz_class, accumulator::Witness>> res = compute_individual_witnesses(witness_multiple, {5, 7, 13});
    classgroup cl;
    for(auto x : res){
        assert(cl.eq(acc.value, cl.exp(x.second.witness.value, x.first)));
    }
}















