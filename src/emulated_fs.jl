using SHA, Random

export FS, get_field, get_query

struct FS
    rng::MersenneTwister 
end 

function FS(seed::Int)
    rng = MersenneTwister(seed)
    return FS(rng)
end

function get_field(fs::FS, ::Type{T}) where T <: BinaryElem
    return rand(fs.rng, T)
end

function get_query(fs::FS, N::Int)
    return rand(fs.rng, 1:N)
end

# mutable struct Transcript
#     state::Vector{UInt8}
# end

# function Transcript()
#     return Transcript(UInt8[])
# end

# TODO: implement this later
# function absorb!(tr::Transcript, data...)
#     for d in data
#         append!(tr.state, reinterpret(UInt8, d))
#     end
# end

# function reseed_rng(tr::Transcript)
#     seed = sha256(tr.state)
#     return MersenneTwister(seed)
# end
