using BinaryFields, MultilinearPoly, Sumcheck, BinaryReedSolomon, MerkleTree

function prover(config::ProverConfig{T, U}, poly::Vector{T}) where {T <: BinaryElem, U <: BinaryElem}
    # initialize fiat shamir emulator: 
    fs = FS(1234)

    # for now just hc this too 
    S = 148

    # initialize proof: 
    proof = LigeritoProof{T, U}()

    println("starting with a first cm")
    # here we commit to 2^24 via matrix of 2^18 * 2^6
    wtns1 = ligero_commit(poly, config.initial_dims[1], config.initial_dims[2], config.initial_reed_solomon)
    println("first cm done")
    cm1 = RecursiveLigeroCommitment(get_root(wtns1.tree))
    # write to the proof
    proof.initial_ligero_cm = cm1

    # then we get first k partial evals at once, before we start the sumcheck
    partial_evals_1 = [get_field(fs, U) for _ in 1:config.initial_k]

    # first time we need to up-cast 
    converted = convert.(U, poly)
    f = MultiLinearPoly(converted)

    # now partially evaluate poly in first k challenges 
    # we don't need to store previous values of f!
    f = partial_eval(f, partial_evals_1) # f is now 2^18

    println("starting with a second cm")
    # now instead of sending it to the verifier let's commit to it again and induce a sumcheck instance
    # this takes 2^18 poly and makes it 2^14 * 2^4 end encodes it
    wtns2 = ligero_commit(f.evals, config.dims[1][1], config.dims[1][2], config.reed_solomon_codes[1])
    cm2 = RecursiveLigeroCommitment(get_root(wtns2.tree))
    push!(proof.recursive_commitments, cm2)
    println("second cm done")
    
    # after committing to this poly we need to induce a sumcheck, so sample random rows and separation challenge 
    rows = size(wtns1.mat, 1)
    queries = sort([get_query(fs, rows) for _ in 1:S])
    alpha = get_field(fs, U)

    
    # TODO! add this to the prover config: 
    sks_vks = eval_sk_at_vks(2^f.n, T)

    println("starting with a merkle opening proof")
    opened_rows = [vec(wtns1.mat[q, :]) for q in queries]
    mtree_proof = MerkleTree.prove(wtns1.tree, queries)
    proof.initial_ligero_proof = RecursiveLigeroProof(opened_rows, mtree_proof)
    println("end merkle opening proof")

    println("starting with a sumcheck induction")
    # finally induce the sumcheck polynomial and enforced sum
    basis_poly, enforced_sum = induce_sumcheck_poly(f.n, sks_vks, opened_rows, partial_evals_1, queries, alpha)
    inner_product = sum(f.evals .* basis_poly)
    @assert inner_product == enforced_sum
    sumcheck_prover = SumcheckProverInstance(f, MultiLinearPoly(basis_poly), enforced_sum)    
    println("end sumcheck induction")
    # now we need to run partial sumcheck for 2^18 -> 2^14, i.e. k[1] rounds of sumcheck
    # then do the gluing that makes sure that 2^18 -> 2^14 step is correct 
    # then run the loop again for 2^14 -> 2^10, i.e. k[2] rounds of sumcheck 

    # then we don't need to send another gluing poly because verifier can check themselves that 2^14 -> 2^10 step is correct
    wtns_prev = wtns2
    cm_prev = cm2
    for i in 1:config.recursive_steps
        println("getting into the loop")
        rs = Vector{U}(undef, config.ks[i])
        for k in 1:config.ks[i]
            ri = get_field(fs, U) 
            fold!(sumcheck_prover, ri)
            rs[k] = ri
        end

        if i == config.recursive_steps
            println("f dim:", sumcheck_prover.f.n)
            println("last step, run ligero")

            # TODO! here we first "send" yr to the verifier before sampling challenges

            rows = size(wtns_prev.mat, 1)
            queries = sort([get_query(fs, rows) for _ in 1:S])

            opened_rows = [vec(wtns_prev.mat[q, :]) for q in queries]
            mtree_proof = MerkleTree.prove(wtns_prev.tree, queries)


            p_final = FinalLigeroProof(sumcheck_prover.f.evals, opened_rows, mtree_proof)
            proof.final_ligero_proof = p_final
            proof.sumcheck_transcript = SumcheckTranscript(sumcheck_prover.transcript)
            return proof
        end 

        wtns_i = ligero_commit(f.evals, config.dims[i][1], config.dims[i][2], config.reed_solomon_codes[i])
        cm_i = RecursiveLigeroCommitment(get_root(wtns_i.tree))
        push!(proof.recursive_commitments, cm_i)

        rows = size(wtns_prev.mat, 1)
        queries = sort([get_query(fs, rows) for _ in 1:S])
        alpha = get_field(fs, U)

        # TODO! add this to the prover config:
        sks_vks = eval_sk_at_vks(2^sumcheck_prover.f.n, U)

        opened_rows = [vec(wtns_prev.mat[q, :]) for q in queries]
        mtree_proof = MerkleTree.prove(wtns_prev.tree, queries)
        p_i = RecursiveLigeroProof(opened_rows, mtree_proof)
        push!(proof.recursive_proofs, p_i)

        basis_poly, enforced_sum = induce_sumcheck_poly(sumcheck_prover.f.n, sks_vks, opened_rows, rs, queries, alpha)
        inner_product = sum(sumcheck_prover.f.evals .* basis_poly)
        @assert inner_product == enforced_sum

        introduce_new!(sumcheck_prover, MultiLinearPoly(basis_poly), enforced_sum)
        beta = get_field(fs, U)
        glue!(sumcheck_prover, beta)

        wtns_prev = wtns_i
        cm_prev = cm_i
    end

    # recursive loop: 
        # for next k_i steps 
        # 0. receive a next evaluation point 
        # 1. run the folding step of the sumcheck 

        # after receiving last challenge eval: 

        # if we are at the last step run normal ligero, 
        # else: 

            # 0. commit to next claimed partial evaluation 
            # 1. receive S queries for the commitment 
            # 2. induce a gluing sumcheck poly based on queries 
            # 3. run that again for the full l recursive steps

    # run the normal ligero + final round of sumcheck 
end

export prover