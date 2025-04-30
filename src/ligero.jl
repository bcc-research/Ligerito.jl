import Base: sizeof
using StatsBase

using BinaryFields, MerkleTree, BinaryReedSolomon
using Base.Threads, ThreadTools

export encode_poly, ligero_commit, verify_ligero

function poly2mat(poly::Vector{T}, m::Int, n::Int, inv_rate::Int) where {T <: BinaryElem}
    m_target = m * inv_rate
    mat = zeros(T, m_target, n)

    nt = Threads.nthreads()
    chunk_size = ceil(Int, n / nt)

    Threads.@sync for t in 1:nt
        Threads.@spawn begin
            start_col = (t-1)*chunk_size + 1
            end_col = min(t*chunk_size, n)

            @inbounds for j in start_col:end_col
                for i in 1:m
                    mat[i, j] = poly[(j-1)*m + i]
                end
            end
        end
    end

    return mat
end

function encode_cols!(poly_mat, rs; parallel=true)
    n = size(poly_mat, 2)

    if parallel
        nt = Threads.nthreads()
        chunk_size = ceil(Int, n / nt)

        Threads.@sync for t in 1:nt
            Threads.@spawn begin
                start_col = (t-1)*chunk_size + 1
                end_col = min(t*chunk_size, n)

                @inbounds for j in start_col:end_col
                    encode_non_systematic!(rs, @view poly_mat[:, j])
                end
            end
        end
    else
        @inbounds for j in 1:n
            encode_non_systematic!(rs, @view poly_mat[:, j])
        end
    end

    return poly_mat
end

function ligero_commit(poly::Vector{T}, m::Int, n::Int, rs::BinaryReedSolomon.ReedSolomonEncoding{T}) where T <: BinaryElem
    poly_mat = poly2mat(poly, m, n, 4)
    encode_cols!(poly_mat, rs)

    leaves = eachrow(poly_mat)
    tree = build_merkle_tree(leaves)

    return RecursiveLigeroWitness(poly_mat, tree)
end

function verify_ligero(queries, opened_rows, yr::Vector{T}, challenges::Vector{U}) where {T, U <: BinaryElem}
    gr = evaluate_lagrange_basis(challenges)
    n = Int(log2(length(yr)))
    sks_vks = eval_sk_at_vks(2^n, U)

    n_threads = Threads.nthreads()
    n_rows = length(opened_rows)
    chunk_size = ceil(Int, n_rows / n_threads)

    Threads.@sync for t in 1:n_threads
        Threads.@spawn begin
            start_idx = (t-1)*chunk_size + 1
            end_idx = min(t*chunk_size, n_rows)

            local_basis = zeros(U, 2^n)
            local_sks_x = Vector{T}(undef, length(sks_vks))  

            @inbounds for i in start_idx:end_idx
                row = opened_rows[i]
                query = queries[i]

                dot = row' * gr

                qf = T(query - 1)
                evaluate_scaled_basis_inplace!(local_sks_x, local_basis, sks_vks, qf, T(1))
                e = yr' * local_basis

                @assert e == dot "Verification failed at index $i: expected $dot, got $e"
            end
        end
    end

end