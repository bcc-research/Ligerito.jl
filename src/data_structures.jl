using BinaryFields, MerkleTree, BinaryReedSolomon

struct ProverConfig{T, U <: BinaryElem}
    recursive_steps::Int
    initial_dims::Tuple{Int, Int}
    dims::Vector{Tuple{Int, Int}} 
    initial_k::Int
    ks::Vector{Int}
    initial_reed_solomon::BinaryReedSolomon.ReedSolomonEncoding{T}
    reed_solomon_codes::Vector{BinaryReedSolomon.ReedSolomonEncoding{U}}
end

struct RecursiveLigeroWitness{T<:BinaryElem}
    mat::Matrix{T}
    tree::CompleteMerkleTree
end

struct RecursiveLigeroCommitment 
    root::MerkleRoot
end 

struct RecursiveLigeroProof{T<:BinaryElem}
    opened_rows::Vector{Vector{T}}
    merkle_proof::BatchedMerkleProof
end

# Note: In case we run just one step of ligero then yr and rows are not always in the same filed, for now we don't deal with this
struct FinalLigeroProof{T<:BinaryElem}
    yr::Vector{T}
    opened_rows::Vector{Vector{T}}
    merkle_proof::BatchedMerkleProof
end

function verify(proof::RecursiveLigeroProof{T}, cm:: RecursiveLigeroCommitment, sorted_queries, depth) where T <: BinaryElem 
    return MerkleTree.verify(cm.root, proof.merkle_proof; depth = depth, leaves = proof.opened_rows, leaf_indices = sorted_queries)
end

struct SumcheckTranscript{T<:BinaryElem}
    tr::Vector{NTuple{3, T}} # at each step prover sends a quadratic polynomial S_i(X)
end

mutable struct LigeritoProof{T, U <: BinaryElem}
    initial_ligero_cm::Union{RecursiveLigeroCommitment, Nothing}
    initial_ligero_proof::Union{RecursiveLigeroProof{T}, Nothing}
    recursive_commitments::Vector{RecursiveLigeroCommitment}
    recursive_proofs::Vector{RecursiveLigeroProof{U}}
    final_ligero_proof::Union{FinalLigeroProof{U}, Nothing}
    sumcheck_transcript::Union{SumcheckTranscript{U}, Nothing}
end

function LigeritoProof{T, U}() where {T, U <: BinaryElem}
    initial_proof = nothing::Union{RecursiveLigeroProof{T}, Nothing} # we need to stabilize T here
    return LigeritoProof{T, U}(
        nothing,
        initial_proof,
        RecursiveLigeroCommitment[],
        Vector{RecursiveLigeroProof{U}}(),  
        nothing,
        nothing
    )
end

export ProverSetup, RecursiveLigeroCommitment, RecursiveLigeroProof, LigeritoProof, SumcheckTranscript, FinalLigeroProof