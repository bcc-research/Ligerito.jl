import Base: sizeof
using StatsBase

using BinaryFields, MerkleTree, BinaryReedSolomon
using Base.Threads, ThreadTools

export encode_poly, ligero_commit


# XXX: Make this allocate way less (a lot of GC pressure)
function encode_cols(poly_mat, rs; parallel=true)
    if parallel
        encoded_columns = tmap(c -> encode_non_systematic!(rs, c), eachcol(poly_mat))
    else
        encoded_columns = map(c -> encode_non_systematic!(rs, c), eachcol(poly_mat))
    end
    return hcat(encoded_columns...)
end

# function encode_cols(poly_mat, rs; parallel=true)
#     n_cols = size(poly_mat, 2)
#     out = Matrix{eltype(poly_mat)}(undef, block_length(rs), n_cols)

#     if parallel
#         Threads.@threads for j in 1:n_cols
#             out[:, j] = encode_non_systematic!(rs, view(poly_mat, :, j))
#         end
#     else
#         for j in 1:n_cols
#             out[:, j] = encode_non_systematic!(rs, view(poly_mat, :, j))
#         end
#     end

#     return out
# end

function ligero_commit(poly::Vector{T}, m::Int, n::Int, rs::BinaryReedSolomon.ReedSolomonEncoding{T}) where T <: BinaryElem
    poly_mat = reshape(poly, m, n)
    mat = encode_cols(poly_mat, rs)

    leaves = eachrow(mat)
    tree = build_merkle_tree(leaves)

    return RecursiveLigeroWitness(mat, tree)
end
