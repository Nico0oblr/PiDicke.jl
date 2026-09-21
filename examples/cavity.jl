using PIDicke

model = DickeModel(8; local_decay = 0.01, local_pump = 0.2)
system = cavity_liouvillian(model; cutoff = 4, g = 0.1, kappa = 1.0)

println("spin dimension: ", length(system.spin))
println("spin-cavity dimension: ", length(system.basis))
println("Liouvillian size: ", size(system.liouvillian.data))
