function _full_dicke_indices(space)
    indices = Dict{Tuple{Int,Int},Int}()
    ranges = block_ranges(space)
    for (sector_index, basis) in enumerate(space.bases)
        S = Int(basis.spinnumber)
        for (local_index, full_index) in enumerate(ranges[sector_index])
            indices[(S, S - (local_index - 1))] = full_index
        end
    end
    return indices
end

function _population_block(liouvillian, space, states)
    indices = _full_dicke_indices(space)
    dimension = length(space)
    block = zeros(ComplexF64, length(states), length(states))
    for (column, source) in enumerate(states)
        operator = zeros(ComplexF64, dimension, dimension)
        source_index = indices[source]
        operator[source_index, source_index] = 1
        derivative = reshape(liouvillian.data * vec(operator), dimension, dimension)
        for (row, target) in enumerate(states)
            target_index = indices[target]
            block[row, column] = derivative[target_index, target_index]
        end
    end
    return block
end

function _coherence_block(liouvillian, space, states)
    indices = _full_dicke_indices(space)
    dimension = length(space)
    block = zeros(ComplexF64, length(states), length(states))
    for (column, (S, M)) in enumerate(states)
        operator = zeros(ComplexF64, dimension, dimension)
        operator[indices[(S, M - 1)], indices[(S, M)]] = 1
        derivative = reshape(liouvillian.data * vec(operator), dimension, dimension)
        for (row, (target_S, target_M)) in enumerate(states)
            block[row, column] = derivative[
                indices[(target_S, target_M - 1)],
                indices[(target_S, target_M)],
            ]
        end
    end
    return block
end

function _expect_after_liouvillian(liouvillian, initial, observable, time)
    dimension = size(initial.data, 1)
    evolved = reshape(
        exp(Matrix(liouvillian.data) * time) * vec(initial.data),
        dimension,
        dimension,
    )
    return real(tr(observable.data * evolved))
end

@testset "cross-representation generator blocks" begin
    model = DickeModel(
        4;
        collective_decay = 0.7,
        collective_pump = 0.11,
        local_decay = 0.13,
        local_pump = 0.19,
        local_dephasing = 0.17,
    )
    full = full_liouvillian(model)

    population_states, _, population_matrix = population_generator(model)
    extracted_population = _population_block(
        full.liouvillian,
        full.space,
        population_states,
    )
    @test isapprox(extracted_population, population_matrix; atol = 1e-12, rtol = 0)

    coherence_states, _, coherence_matrix = first_order_coherence_generator(model)
    extracted_coherence = _coherence_block(
        full.liouvillian,
        full.space,
        coherence_states,
    )
    @test isapprox(extracted_coherence, coherence_matrix; atol = 1e-12, rtol = 0)
end

@testset "bad-cavity convention" begin
    N = 1
    g = 1.0
    kappa = 40.0
    collective_rate = 4g^2 / kappa
    scaled_times = range(0.25, 4.0, length = 16)
    times = scaled_times ./ collective_rate

    cavity_model = DickeModel(N)
    cavity = cavity_liouvillian(cavity_model; cutoff = 1, g, kappa)
    cavity_initial = dm(fully_excited_state(cavity_model) ⊗ fockstate(cavity.cavity, 0))
    cavity_inversion = [
        _expect_after_liouvillian(cavity.liouvillian, cavity_initial, cavity.Sz, time)
        for time in times
    ]

    eliminated_model = DickeModel(N; collective_decay = collective_rate)
    eliminated = full_liouvillian(eliminated_model)
    eliminated_initial = dm(fully_excited_state(eliminated_model))
    eliminated_Sz = collective_operator(eliminated_model, :z)
    eliminated_inversion = [
        _expect_after_liouvillian(
            eliminated.liouvillian,
            eliminated_initial,
            eliminated_Sz,
            time,
        )
        for time in times
    ]

    @test maximum(abs, cavity_inversion - eliminated_inversion) < 0.01
end
