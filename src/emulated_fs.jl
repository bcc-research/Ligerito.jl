using SHA, Random, Sumcheck

export FS, get_field, get_query, absorb!

# struct FS
#     rng::MersenneTwister 
# end 

# function FS(seed::Int)
#     rng = MersenneTwister(seed)
#     return FS(rng)
# end

# function get_field(fs::FS, ::Type{T}) where T <: BinaryElem
#     return rand(fs.rng, T)
# end

# function get_query(fs::FS, N::Int)
#     return rand(fs.rng, 1:N)
# end

# the most basic implementation of Fiat-Shamir
mutable struct FS
    ctx::SHA256_CTX
    counter::UInt32
end

function FS(seed::Int)
    ctx = SHA.SHA256_CTX()
    SHA.update!(ctx, reinterpret(UInt8, [seed]))
    return FS(ctx, 0)
end

function absorb!(fs::FS, root::MerkleRoot)
    SHA.update!(fs.ctx, root.root)
end

function absorb!(fs::FS, elems::Vector{T}) where T <: BinaryElem
    SHA.update!(fs.ctx, reinterpret(UInt8, elems))
end

function absorb!(fs::FS, s::QuadraticEvals)
    SHA.update!(fs.ctx, reinterpret(UInt8, [s.e0, s.e1, s.e2]))
end

function squeeze_rng!(fs::FS)::MersenneTwister
    SHA.update!(fs.ctx, reinterpret(UInt8, [fs.counter]))
    fs.counter += 1
    digest = SHA.digest!(deepcopy(fs.ctx))
    return MersenneTwister(reinterpret(UInt32, digest))
end

function get_field(fs::FS, ::Type{T}) where T <: BinaryElem
    rng = squeeze_rng!(fs)
    return rand(rng, T)
end

function get_query(fs::FS, N::Int)
    rng = squeeze_rng!(fs)
    return rand(rng, 1:N)
end

function get_distinct_queries(fs::FS, N::Int, S::Int)
    queries = Int[]
    seen = Set{Int}()

    while length(queries) < S
        q = get_query(fs, N)
        if q ∉ seen
            push!(queries, q)
            push!(seen, q)
        end
    end

    return sort(queries)
end
