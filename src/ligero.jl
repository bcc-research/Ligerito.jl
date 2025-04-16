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
    @time mat = encode_cols(poly_mat, rs)

    leaves = eachrow(mat)
    @time tree = build_merkle_tree(leaves)

    return RecursiveLigeroWitness(mat, tree)
end
