using BinaryFields

export verifier

function verifier(proof::FinalizedLigeritoProof{T, U}) where {T <: BinaryElem, U <: BinaryElem}
    fs = FS(1234)
    S = 148

    absorb!(fs, proof.initial_ligero_cm.root)
    partial_evals_1 = [get_field(fs, U) for _ in 1:6]

    absorb!(fs, proof.recursive_commitments[1].root)

    queries = get_distinct_queries(fs, 2^20, S)
    # given queries we can now verify the merkle tree proof 
    res = MerkleTree.verify(proof.initial_ligero_cm.root, proof.initial_ligero_proof.merkle_proof; depth = 20, leaves = proof.initial_ligero_proof.opened_rows, leaf_indices = queries)
    @assert res == true

    alpha = get_field(fs, U)

    # TODO! add this to the verifier config: 
    sks_vks = eval_sk_at_vks(2^18, T)
    basis_poly, enforced_sum = induce_sumcheck_poly_parallel(18, sks_vks, proof.initial_ligero_proof.opened_rows, partial_evals_1, queries, alpha) 

    sumcheck_verifier, g1 = SumcheckVerifierInstance(MultiLinearPoly(basis_poly), enforced_sum, proof.sumcheck_transcript.tr)
    absorb!(fs, g1)
    for i in 1:2
        rs = Vector{U}(undef, 4)
        for k in 1:4
            ri = get_field(fs, U) 
            si = fold!(sumcheck_verifier, ri)
            absorb!(fs, si)
            rs[k] = ri
        end

        root = proof.recursive_commitments[i].root

        if i == 2
            absorb!(fs, proof.final_ligero_proof.yr)

            queries = get_distinct_queries(fs, 2^16, S)
            res = MerkleTree.verify(root, proof.final_ligero_proof.merkle_proof; depth = 16, leaves = proof.final_ligero_proof.opened_rows, leaf_indices = queries)
            @assert res == true

            # take last random element
            final_r = get_field(fs, U)
            f = MultiLinearPoly(proof.final_ligero_proof.yr)
            f_eval = partial_eval(f, [final_r]).evals

            ok = Sumcheck.verify_partial(sumcheck_verifier, final_r, f_eval)
            return ok
        end 

        liger_proof = proof.recursive_proofs[i]
        absorb!(fs, proof.recursive_commitments[i + 1].root)
        queries = get_distinct_queries(fs, 2^16, S)
        res = MerkleTree.verify(root, liger_proof.merkle_proof; depth = 16, leaves = liger_proof.opened_rows, leaf_indices = queries)
        @assert res == true


        alpha = get_field(fs, U)

        # TODO! add this to the verifier config: 
        sks_vks = eval_sk_at_vks(2^14, U)
        basis_poly, enforced_sum = induce_sumcheck_poly_parallel(14, sks_vks, liger_proof.opened_rows, rs, queries, alpha)

        gl_i = introduce_new!(sumcheck_verifier, MultiLinearPoly(basis_poly), enforced_sum)
        absorb!(fs, gl_i)

        beta = get_field(fs, U)
        glue!(sumcheck_verifier, beta)
    end
end