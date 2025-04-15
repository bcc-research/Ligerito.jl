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
    opened_rows::Vector{Vector{T}} # each row is a vector of binary elements
    merkle_proof::BatchedMerkleProof
end

function verify(proof::RecursiveLigeroProof{T}, cm:: RecursiveLigeroCommitment, sorted_queries, depth) where T <: BinaryElem 
    return MerkleTree.verify(cm.root, proof.merkle_proof; depth = depth, leaves = proof.opened_rows, leaf_indices = sorted_queries)
end

struct SumcheckTranscript{T<:BinaryElem}
    tr::Vector{NTuple{3, T}} # at each step prover sends a quadratic polynomial S_i(X)
end

struct LigeritoProof{T <: BinaryElem}
    recursive_commitments::Vector{RecursiveLigeroCommitment}    
    recursive_proofs::Vector{RecursiveLigeroProof{T}}
    sumcheck_transcript::SumcheckTranscript{T}
end

export ProverSetup, RecursiveLigeroCommitment, RecursiveLigeroProof, LigeritoProof, SumcheckTranscript 