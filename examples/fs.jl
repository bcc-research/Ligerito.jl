using BinaryFields, SHA

# fs = FS(1234)
x = rand(BinaryElem32, 10)
@show sizeof(x)

y = [copy(x), copy(x), copy(x), copy(x)]
@show sizeof(y)

# absorb!(fs, x)
# s = get_field(fs, BinaryElem128)
# @show s
# s = get_field(fs, BinaryElem128)
# @show s

# x = get_field(fs, BinaryElem128)
# print(x)

# # reinterpret(UInt8, x)

# res = sha256(reinterpret(UInt8, x))
# @show res

# z::Int = 1234 
# @show reinterpret(UInt8, z) 