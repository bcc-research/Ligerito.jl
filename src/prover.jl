using BinaryFields, MultilinearPoly, Sumcheck, BinaryReedSolomon, MerkleTree
using StatsBase

function convert_parallel!(converted::Vector{U}, poly::Vector{T}) where {U, T}
    n = length(poly)
    nt = Threads.nthreads()
    chunk_size = ceil(Int, n / nt)

    Threads.@sync for t in 1:nt
        Threads.@spawn begin
            start_idx = (t - 1) * chunk_size + 1
            end_idx = min(t * chunk_size, n)
            @inbounds for i in start_idx:end_idx
                converted[i] = convert(U, poly[i])
            end
        end
    end
end

function convert_x(x::Vector{T}) where T <: BinaryElem
    converted = Vector{BinaryElem128}(undef, length(x))
    Threads.@threads for i in 1:length(x)
        converted[i] = convert(BinaryElem128, x[i])
    end
    return converted    
end

function prover(config::ProverConfig{T, U}, poly::Vector{T}) where {T <: BinaryElem, U <: BinaryElem}
    # initialize fiat shamir: 
    fs = FS(1234)

    # for now just hc this too 
    S = 148

    # initialize proof: 
    proof = LigeritoProof{T, U}()

    # here we commit to 2^24 via matrix of 2^18 * 2^6
    wtns1 = ligero_commit(poly, config.initial_dims[1], config.initial_dims[2], config.initial_reed_solomon)
    cm1 = RecursiveLigeroCommitment(get_root(wtns1.tree))
    proof.initial_ligero_cm = cm1
    absorb!(fs, cm1.root)
    # then we get first k partial evals at once, before we start the sumcheck
    partial_evals_1 = [get_field(fs, U) for _ in 1:config.initial_k]

    # first time we need to up-cast 
    # @time begin
    #     converted = Vector{U}(undef, length(poly))
    #     Threads.@threads for i in 1:length(poly)
    #         converted[i] = convert(U, poly[i])
    #     end
    # end 
    converted = convert_x(poly)
    # converted = Vector{U}(undef, length(poly))
    # @time convert_parallel!(converted, poly)
    f = MultiLinearPoly(converted)

    # now partially evaluate poly in first k challenges 
    # we don't need to store previous values of f!
    f = partial_eval(f, partial_evals_1) # f is now 2^18

    # now instead of sending it to the verifier let's commit to it again and induce a sumcheck instance
    # this takes 2^18 poly and makes it 2^14 * 2^4 end encodes it
    wtns2 = ligero_commit(f.evals, config.dims[1][1], config.dims[1][2], config.reed_solomon_codes[1])
    cm2 = RecursiveLigeroCommitment(get_root(wtns2.tree))
    push!(proof.recursive_commitments, cm2)
    absorb!(fs, cm2.root)
    
    # after committing to this poly we need to induce a sumcheck, so sample random rows and separation challenge 
    rows = size(wtns1.mat, 1)
    queries = get_distinct_queries(fs, rows, S)
    alpha = get_field(fs, U)

    # TODO! add this to the prover config: 
    sks_vks = eval_sk_at_vks(2^f.n, T)

    opened_rows = [vec(wtns1.mat[q, :]) for q in queries]
    mtree_proof = MerkleTree.prove(wtns1.tree, queries)
    proof.initial_ligero_proof = RecursiveLigeroProof(opened_rows, mtree_proof)

    # finally induce the sumcheck polynomial and enforced sum
    basis_poly, enforced_sum = induce_sumcheck_poly_parallel(f.n, sks_vks, opened_rows, partial_evals_1, queries, alpha) 
    inner_product = sum(f.evals .* basis_poly)
    @assert inner_product == enforced_sum
    sumcheck_prover, s1 = SumcheckProverInstance(f, MultiLinearPoly(basis_poly), enforced_sum)   
    absorb!(fs, s1) 
    # now we need to run partial sumcheck for 2^18 -> 2^14, i.e. k[1] rounds of sumcheck
    # then do the gluing that makes sure that 2^18 -> 2^14 step is correct 
    # then run the loop again for 2^14 -> 2^10, i.e. k[2] rounds of sumcheck 
    # then we don't need to send another gluing poly because verifier can check themselves that 2^14 -> 2^10 step is correct

    wtns_prev = wtns2
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
            queries = get_distinct_queries(fs, rows, S)
            opened_rows = [vec(wtns_prev.mat[q, :]) for q in queries]
            mtree_proof = MerkleTree.prove(wtns_prev.tree, queries)


            p_final = FinalLigeroProof(sumcheck_prover.f.evals, opened_rows, mtree_proof)
            proof.final_ligero_proof = p_final
            proof.sumcheck_transcript = SumcheckTranscript(sumcheck_prover.transcript)
            return finalize(proof)
        end 

        wtns_i = ligero_commit(f.evals, config.dims[i][1], config.dims[i][2], config.reed_solomon_codes[i])
        cm_i = RecursiveLigeroCommitment(get_root(wtns_i.tree))
        push!(proof.recursive_commitments, cm_i)
        absorb!(fs, cm_i.root)

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
        inner_product = sum(sumcheck_prover.f.evals .* basis_poly)
        @assert inner_product == enforced_sum

        gl_i = introduce_new!(sumcheck_prover, MultiLinearPoly(basis_poly), enforced_sum)
        absorb!(fs, gl_i)

        beta = get_field(fs, U)
        glue!(sumcheck_prover, beta)

        wtns_prev = wtns_i
    end
end

export prover