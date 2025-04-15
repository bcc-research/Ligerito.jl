using BinaryFields

config = hardcoded_config(BinaryElem32, BinaryElem128)
poly = rand(BinaryElem32, 2^24)

prover(config, poly)
