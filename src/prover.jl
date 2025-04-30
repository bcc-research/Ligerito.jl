using BinaryFields, MultilinearPoly, Sumcheck, BinaryReedSolomon, MerkleTree
using StatsBase

# S is hardcoded to 148 for now 
# Rate is hardcoded to 1/4 for now
function prover(config::ProverConfig{T, U}, poly::Vector{T}) where {T <: BinaryElem, U <: BinaryElem}
    # initialize fiat shamir: 
    # TODO! hash public data into it 
    fs = FS(1234)

    # for now just hc this too 
    S = 148

    proof = LigeritoProof{T, U}()

    wtns_0 = ligero_commit(poly, config.initial_dims[1], config.initial_dims[2], config.initial_reed_solomon)
    cm_0 = RecursiveLigeroCommitment(get_root(wtns_0.tree))
    proof.initial_ligero_cm = cm_0
    absorb!(fs, cm_0.root)

    partial_evals_0 = [get_field(fs, U) for _ in 1:config.initial_k]

    f = MultiLinearPoly(poly)
    f = partial_eval(f, partial_evals_0)

    wtns_1 = ligero_commit(f.evals, config.dims[1][1], config.dims[1][2], config.reed_solomon_codes[1])
    cm_1 = RecursiveLigeroCommitment(get_root(wtns_1.tree))
    push!(proof.recursive_commitments, cm_1)
    absorb!(fs, cm_1.root)
    
    rows = size(wtns_0.mat, 1)
    queries = get_distinct_queries(fs, rows, S)
    alpha = get_field(fs, U)

    # TODO! add this to the prover config: 
    sks_vks = eval_sk_at_vks(2^f.n, T)

    opened_rows = [vec(wtns_0.mat[q, :]) for q in queries]
    mtree_proof = MerkleTree.prove(wtns_0.tree, queries)
    proof.initial_ligero_proof = RecursiveLigeroProof(opened_rows, mtree_proof)

    basis_poly, enforced_sum = induce_sumcheck_poly_parallel(f.n, sks_vks, opened_rows, partial_evals_0, queries, alpha) 
    sumcheck_prover, s1 = SumcheckProverInstance(f, MultiLinearPoly(basis_poly), enforced_sum)   
    absorb!(fs, s1) 

    wtns_prev = wtns_1
    for i in 1:config.recursive_steps
        rs = Vector{U}(undef, config.ks[i])
        for k in 1:config.ks[i]
            ri = get_field(fs, U) 
            si = fold!(sumcheck_prover, ri)
            absorb!(fs, si)
            rs[k] = ri
        end

        if i == config.recursive_steps
            absorb!(fs, sumcheck_prover.f.evals)

            rows = size(wtns_prev.mat, 1)
            # @show rows
            queries = get_distinct_queries(fs, rows, S)
            opened_rows = [vec(wtns_prev.mat[q, :]) for q in queries]
            mtree_proof = MerkleTree.prove(wtns_prev.tree, queries)

            p_final = FinalLigeroProof(sumcheck_prover.f.evals, opened_rows, mtree_proof)
            proof.final_ligero_proof = p_final
            proof.sumcheck_transcript = SumcheckTranscript(sumcheck_prover.transcript)
            return finalize(proof)
        end 

        wtns_next = ligero_commit(sumcheck_prover.f.evals, config.dims[i + 1][1], config.dims[i + 1][2], config.reed_solomon_codes[i + 1])
        cm_next = RecursiveLigeroCommitment(get_root(wtns_next.tree))
        push!(proof.recursive_commitments, cm_next)
        absorb!(fs, cm_next.root)

        rows = size(wtns_prev.mat, 1)
        queries = get_distinct_queries(fs, rows, S)
        alpha = get_field(fs, U)

        # TODO! add this to the prover config:
        sks_vks = eval_sk_at_vks(2^sumcheck_prover.f.n, U)

        opened_rows = [vec(wtns_prev.mat[q, :]) for q in queries]
        mtree_proof = MerkleTree.prove(wtns_prev.tree, queries)
        p_i = RecursiveLigeroProof(opened_rows, mtree_proof)
        push!(proof.recursive_proofs, p_i)

        basis_poly, enforced_sum = induce_sumcheck_poly_parallel(sumcheck_prover.f.n, sks_vks, opened_rows, rs, queries, alpha)

        gl_i = introduce_new!(sumcheck_prover, MultiLinearPoly(basis_poly), enforced_sum)
        absorb!(fs, gl_i)

        beta = get_field(fs, U)
        glue!(sumcheck_prover, beta)

        wtns_prev = wtns_next
    end
end

export prover