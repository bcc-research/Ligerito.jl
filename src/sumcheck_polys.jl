using BinaryFields

function precompute_alpha_powers(α::T, n::Int) where {T}
    α_pows = Vector{T}(undef, n)
    α_pows[1] = one(T)

    for i in 2:n
        α_pows[i] = α_pows[i-1] * α
    end

    return α_pows
end

function induce_sumcheck_poly_parallel(n::Int, sks_vks::Vector{T}, opened_rows::Vector{Vector{T}}, v_challenges::Vector{U}, sorted_queries::Vector{Int}, α::U) where {U <: BinaryElem, T <: BinaryElem}
    gr = evaluate_lagrange_basis(v_challenges)
    @assert all(length(row) == length(gr) for row in opened_rows)
    @assert length(opened_rows) == length(sorted_queries)

    n_threads = Threads.nthreads()
    partial_basis = [zeros(U, 2^n) for _ in 1:n_threads]
    partial_sums  = zeros(U, n_threads)

    n_rows = length(opened_rows)
    chunk_size = ceil(Int, n_rows / n_threads)

    alpha_pows = precompute_alpha_powers(α, n_rows)

    Threads.@sync for t in 1:n_threads
        Threads.@spawn begin
            local_basis = zeros(U, 2^n)
            local_sks_x = Vector{T}(undef, length(sks_vks))   

            start_idx = (t-1)*chunk_size + 1
            end_idx = min(t*chunk_size, n_rows)

            @inbounds for i in start_idx:end_idx
                row = opened_rows[i]
                query = sorted_queries[i]

                dot = row' * gr

                α_pow = alpha_pows[i]

                #NOTE: This is not needed on a prover side, but it's here for sanity checks
                partial_sums[t] += dot * α_pow

                qf = T(query - 1)
                evaluate_scaled_basis_inplace!(local_sks_x, local_basis, sks_vks, qf, α_pow)
                @. partial_basis[t] += local_basis
            end
        end
    end

    basis_poly = sum(partial_basis)
    enforced_sum = sum(partial_sums)

    return basis_poly, enforced_sum
end

export induce_sumcheck_poly, evaluate_lagrange_basis, induce_sumcheck_poly_parallel