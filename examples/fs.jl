using BinaryFields

fs = FS(1234)
x = get_field(fs, BinaryElem128)
print(x)