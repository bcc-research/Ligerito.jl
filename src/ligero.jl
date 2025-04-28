import Base: sizeof
using StatsBase

using BinaryFields, MerkleTree, BinaryReedSolomon
using Base.Threads, ThreadTools

export encode_poly, ligero_commit

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
    tree = build_merkle_tree_fast(leaves)

    return RecursiveLigeroWitness(poly_mat, tree)
end
