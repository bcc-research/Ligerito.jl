using BinaryFields

export verifier

#TODO: Optimize verifier with time and allocations!

function verifier(config::VerifierConfig, proof::FinalizedLigeritoProof{T, U}) where {T <: BinaryElem, U <: BinaryElem}
    fs = FS(1234)
    S = 148
    log_inv_rate = 2

    absorb!(fs, proof.initial_ligero_cm.root)
    partial_evals_0 = [get_field(fs, U) for _ in 1:config.initial_k]

    absorb!(fs, proof.recursive_commitments[1].root)

    depth = config.initial_dim + log_inv_rate
    queries = get_distinct_queries(fs, 2^depth, S)
    res = MerkleTree.verify(proof.initial_ligero_cm.root, proof.initial_ligero_proof.merkle_proof; depth, leaves = proof.initial_ligero_proof.opened_rows, leaf_indices = queries)
    @assert res == true

    alpha = get_field(fs, U)

    # TODO! add this to the verifier config: 
    sks_vks = eval_sk_at_vks(2^config.initial_dim, T)
    basis_poly, enforced_sum = induce_sumcheck_poly_parallel(config.initial_dim, sks_vks, proof.initial_ligero_proof.opened_rows, partial_evals_0, queries, alpha) 

    sumcheck_verifier, g1 = SumcheckVerifierInstance(MultiLinearPoly(basis_poly), enforced_sum, proof.sumcheck_transcript.tr)
    absorb!(fs, g1)
    for i in 1:config.recursive_steps
        rs = Vector{U}(undef, config.ks[i])
        for k in 1:config.ks[i]
            ri = get_field(fs, U) 
            si = fold!(sumcheck_verifier, ri)
            absorb!(fs, si)
            rs[k] = ri
        end

        root = proof.recursive_commitments[i].root
        if i == config.recursive_steps
            absorb!(fs, proof.final_ligero_proof.yr)

            depth = config.log_dims[i] + log_inv_rate
            queries = get_distinct_queries(fs, 2^depth, S)
            res = MerkleTree.verify(root, proof.final_ligero_proof.merkle_proof; depth, leaves = proof.final_ligero_proof.opened_rows, leaf_indices = queries)
            @assert res == true

            verify_ligero(queries, proof.final_ligero_proof.opened_rows, proof.final_ligero_proof.yr, rs)

            # take last random element
            final_r = get_field(fs, U)
            f = MultiLinearPoly(proof.final_ligero_proof.yr)
            f_eval = partial_eval(f, [final_r]).evals

            ok = Sumcheck.verify_partial(sumcheck_verifier, final_r, f_eval)
            return ok
        end 

        absorb!(fs, proof.recursive_commitments[i + 1].root)

        depth = config.log_dims[i] + log_inv_rate
        liger_proof = proof.recursive_proofs[i]
        queries = get_distinct_queries(fs, 2^depth, S)
        res = MerkleTree.verify(root, liger_proof.merkle_proof; depth, leaves = liger_proof.opened_rows, leaf_indices = queries)
        @assert res == true

        alpha = get_field(fs, U)

        # TODO! add this to the verifier config: 
        sks_vks = eval_sk_at_vks(2^config.log_dims[i], U)
        basis_poly, enforced_sum = induce_sumcheck_poly_parallel(config.log_dims[i], sks_vks, liger_proof.opened_rows, rs, queries, alpha)

        gl_i = introduce_new!(sumcheck_verifier, MultiLinearPoly(basis_poly), enforced_sum)
        absorb!(fs, gl_i)

        beta = get_field(fs, U)
        glue!(sumcheck_verifier, beta)
    end
end