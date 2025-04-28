using BinaryFields, MerkleTree, BinaryReedSolomon
import Base: sizeof

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
sizeof(x::RecursiveLigeroCommitment) = sizeof(x.root)   

struct RecursiveLigeroProof{T<:BinaryElem}
    opened_rows::Vector{Vector{T}}
    merkle_proof::BatchedMerkleProof
end
sizeof(x::RecursiveLigeroProof) = sizeof(x.opened_rows) + sizeof(x.merkle_proof)

# Note: In case we run just one step of ligero then yr and rows are not always in the same filed, for now we don't deal with this
struct FinalLigeroProof{T<:BinaryElem}
    yr::Vector{T}
    opened_rows::Vector{Vector{T}}
    merkle_proof::BatchedMerkleProof
end
sizeof(x::FinalLigeroProof) = sizeof(x.yr) + sizeof(x.opened_rows) + sizeof(x.merkle_proof)

function verify(proof::RecursiveLigeroProof{T}, cm:: RecursiveLigeroCommitment, sorted_queries, depth) where T <: BinaryElem 
    return MerkleTree.verify(cm.root, proof.merkle_proof; depth = depth, leaves = proof.opened_rows, leaf_indices = sorted_queries)
end

struct SumcheckTranscript{T<:BinaryElem}
    tr::Vector{NTuple{3, T}} # at each step prover sends a quadratic polynomial S_i(X)
end
Base.sizeof(_::NTuple{3, T}) where {T<:BinaryElem} = 3 * sizeof(T)
Base.sizeof(x::SumcheckTranscript) = sum(sizeof, x.tr)

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

struct FinalizedLigeritoProof{T, U <: BinaryElem}
    initial_ligero_cm::RecursiveLigeroCommitment
    initial_ligero_proof::RecursiveLigeroProof{T}
    recursive_commitments::Vector{RecursiveLigeroCommitment}
    recursive_proofs::Vector{RecursiveLigeroProof{U}}
    final_ligero_proof::FinalLigeroProof{U}
    sumcheck_transcript::SumcheckTranscript{U}
end

function finalize(p::LigeritoProof{T, U}) where {T, U <: BinaryElem}
    @assert p.initial_ligero_cm !== nothing "Missing initial_ligero_cm"
    @assert p.initial_ligero_proof !== nothing "Missing initial_ligero_proof"
    @assert p.final_ligero_proof !== nothing "Missing final_ligero_proof"
    @assert p.sumcheck_transcript !== nothing "Missing sumcheck_transcript"

    return FinalizedLigeritoProof{T, U}(
        p.initial_ligero_cm,
        p.initial_ligero_proof,
        p.recursive_commitments,
        p.recursive_proofs,
        p.final_ligero_proof,
        p.sumcheck_transcript,
    )
end

Base.sizeof(p::FinalizedLigeritoProof{T, U}) where {T, U <: BinaryElem} =
    sizeof(p.initial_ligero_cm) +
    sizeof(p.initial_ligero_proof) +
    sum(sizeof.(p.recursive_commitments)) +
    sum(sizeof.(p.recursive_proofs)) +
    sizeof(p.final_ligero_proof) +
    sizeof(p.sumcheck_transcript)

export ProverSetup, RecursiveLigeroCommitment, RecursiveLigeroProof, LigeritoProof, SumcheckTranscript, FinalLigeroProof
export finalize, FinalizedLigeritoProof