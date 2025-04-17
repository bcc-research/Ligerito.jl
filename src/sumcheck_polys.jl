using BinaryFields

function evaluate_lagrange_basis(rs::Vector{T}) where T <: BinaryElem
    one_elem = one(T)
    current_layer = [one_elem + rs[1], rs[1]]
    len = 2
    for i in 2:length(rs)
        next_layer_size = 2 * len
        next_layer = Vector{T}(undef, next_layer_size)

        ri_p_one = one_elem + rs[i]
        for j in 1:len
            next_layer[2*j - 1] = current_layer[j] * ri_p_one
            next_layer[2*j]   = current_layer[j] * rs[i]
        end

        current_layer = next_layer
        len *= 2
    end

    return current_layer
end

function induce_sumcheck_poly(n::Int, sks_vks::Vector{T}, opened_rows::Vector{Vector{T}}, v_challenges::Vector{U}, sorted_queries::Vector{Int}, α::U) where {U <: BinaryElem, T <: BinaryElem}
    gr = evaluate_lagrange_basis(v_challenges)
    @assert all(length(row) == length(gr) for row in opened_rows)
    @assert length(opened_rows) == length(sorted_queries)


    basis_poly = zeros(U, 2^n)
    enforced_sum = zero(U)
    α_pow = one(U)

    for (row, query) in zip(opened_rows, sorted_queries)
        dot = sum(row .* gr)
        enforced_sum += dot * α_pow

        qf = T(query - 1)
        basis_q_evals = evaluate_basis(2^n, sks_vks, qf)
        basis_poly .+= α_pow .* basis_q_evals

        α_pow *= α
    end

    return (basis_poly, enforced_sum)
end 

function induce_sumcheck_poly_parallel(n::Int, sks_vks::Vector{T}, opened_rows::Vector{Vector{T}}, v_challenges::Vector{U}, sorted_queries::Vector{Int}, α::U) where {U <: BinaryElem, T <: BinaryElem}
    gr = evaluate_lagrange_basis(v_challenges)
    @assert all(length(row) == length(gr) for row in opened_rows)
    @assert length(opened_rows) == length(sorted_queries)

    n_threads = Threads.nthreads()
    partial_basis = [zeros(U, 2^n) for _ in 1:n_threads]
    partial_sums  = [zero(U) for _ in 1:n_threads]

    Threads.@threads for i in 1:length(opened_rows)
        tid = Threads.threadid()
        row = opened_rows[i]
        query = sorted_queries[i]

        dot = sum(row .* gr)

        α_pow = α^(i - 1)  # Later we cna first run all powers

        partial_sums[tid] += dot * α_pow

        qf = T(query - 1)
        basis_q_evals = evaluate_basis(2^n, sks_vks, qf)

        @. partial_basis[tid] += α_pow * basis_q_evals
    end

    basis_poly = sum(partial_basis)
    enforced_sum = sum(partial_sums)

    return basis_poly, enforced_sum
end


export induce_sumcheck_poly, evaluate_lagrange_basis, induce_sumcheck_poly_parallel