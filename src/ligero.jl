import Base: sizeof

using BinaryFields, MerkleTree, BinaryReedSolomon

export encode_poly, ligero_commit

next_pow_two(x) = 1 << (1 + floor(Int, log2(x - 1)))

# n_queries(; log2_prob=100, inv_rate=4) = ceil(Int, -log2_prob/log2((1+1/inv_rate)/2))

# Assumes `N` is a power of 2
# function opt_dims(N, prop::LigeroProofProperties)
#     S = n_queries(; prop.log2_prob, prop.inv_rate)
#     n = next_pow_two(sqrt(prop.ext_field_size*N/(prop.base_field_size*S)))

#     return div(N, n), n
# end

# struct LigeroCommitment{T}
#     mat::Matrix{T}
#     tree::CompleteMerkleTree
#     rs::BinaryReedSolomon.ReedSolomonEncoding{T}
# end

# matrix(c::LigeroCommitment) = c.mat

# struct LigeroVerifierCommitment
#     root::MerkleRoot
# end

# sizeof(x::LigeroVerifierCommitment) = sizeof(x.root)

# verifier_commitment(com) = LigeroVerifierCommitment(get_root(com.tree))

# XXX: Make this allocate way less (a lot of GC pressure)
function encode_cols(poly_mat, rs; parallel=true)
    # if parallel
    #     encoded_columns = tmap(c -> encode(rs, c), eachcol(poly_mat))
    # else
    encoded_columns = map(c -> encode_non_systematic!(rs, c), eachcol(poly_mat))
    # end
    return hcat(encoded_columns...)
end

function encode_poly(poly, properties=nothing)
    T_poly = eltype(poly)

    if isnothing(properties)
        properties = LigeroProofProperties(100, 4, bitsize(T_poly), 128)
    end

    m, n = opt_dims(length(poly), properties)
    rs = reed_solomon(T_poly, m, m*properties.inv_rate)

    poly_mat = reshape(poly, m, n)
    mat = encode_cols(poly_mat, rs)

    mat
end

function ligero_commit(poly::Vector{T}, m::Int, n::Int, rs::BinaryReedSolomon.ReedSolomonEncoding{T}) where T <: BinaryElem
    poly_mat = reshape(poly, m, n)
    mat = encode_cols(poly_mat, rs)

    leaves = eachrow(mat)
    tree = build_merkle_tree(leaves)

    return RecursiveLigeroWitness(mat, tree)
end

# function commit(poly; properties=nothing, verbose=false, parallel=true)
#     @assert length(poly) > 0 "Polynomial must have at least one element"
#     @assert is_pow_two(length(poly)) "Polynomial length must be a power of 2"

#     T_poly = eltype(poly)

#     if isnothing(properties)
#         properties = LigeroProofProperties(100, 4, bitsize(T_poly), 128)
#     end

#     m, n = opt_dims(length(poly), properties)
#     min_field_size = round(Int, log2(m))+round(Int, log2(properties.inv_rate))
#     if min_field_size > properties.base_field_size
#         error("Base field size $(properties.base_field_size) is too small for polynomial of length $(length(poly)) with inv_rate $(properties.inv_rate). Minimum size is $min_field_size.")
#     end
#     rs = reed_solomon(T_poly, m, m*properties.inv_rate)

#     poly_mat = reshape(poly, m, n)
#     if verbose
#         @info "Encoding"
#         @time mat = encode_cols(poly_mat, rs; parallel)
#     else
#         mat = encode_cols(poly_mat, rs)
#     end

#     leaves = eachrow(mat)
    
#     if verbose
#         @info "Building Merkle tree"
#         @time tree = build_merkle_tree(leaves)
#     else
#         tree = build_merkle_tree(leaves)
#     end

#     return LigeroCommitment(mat, tree, rs)
# end
