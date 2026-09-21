using LinearAlgebra
using PIDicke
using QuantumOptics: dm
using Random
using Test

@testset "PIDicke" begin
    @testset "model and bases" begin
        model = DickeModel(4; collective_decay = 1.0)
        @test jmax(model) == 2
        @test length(dicke_triangle(model)) == 9
        @test length(dicke_subspaces(4)) == 3
        @test length(dicke_space(model)) == 9
        @test length(dicke_space(3)) == 6
        @test_throws ArgumentError dicke_triangle(3)
        @test representation_basis(model, PopulationRepresentation()) == dicke_triangle(model)
        @test length(representation_basis(model, FullDensityRepresentation())) == 9
    end

    @testset "transition amplitudes" begin
        @test A_JM_minus(2, 2) == 2
        @test A_JM_plus(2, -2) == 2
        @test P_JM_minus_0(0, 0, 4) == 0
        @test P_JM_z_minus(0, 0, 4) == 0
        @test d_N_J(10, 3) == 1
    end

    @testset "populations and trajectories" begin
        model = DickeModel(
            4;
            collective_decay = 1.0,
            local_decay = 0.1,
            local_pump = 0.2,
            local_dephasing = 0.05,
        )
        states, indices, generator = population_generator(model)
        @test size(generator) == (9, 9)
        @test maximum(abs, vec(sum(generator; dims = 1))) < 1e-12

        initial = population_initial_state(indices, (2, 2))
        @test sum(initial) == 1
        steady = population_steady_state(generator)
        @test isapprox(sum(steady), 1; atol = 1e-12)
        @test norm(generator * steady) < 1e-10

        times, trajectory = run_population_trajectory(
            model;
            tmax = 0.2,
            rng = Xoshiro(1234),
        )
        @test issorted(times)
        @test all(state -> begin
            S, M = state
            0 <= S <= 2 && abs(M) <= S
        end, trajectory)
    end

    @testset "first-order coherence band" begin
        model = DickeModel(4; collective_decay = 1.0, local_pump = 0.2)
        states, indices, generator = first_order_coherence_generator(model)
        @test length(states) == 6
        @test size(generator) == (6, 6)
        @test all(state -> begin
            S, M = state
            abs(M) <= S && abs(M - 1) <= S
        end, states)

        population_states, _, _ = population_generator(model)
        populations = population_steady_state(model)
        initial = coherence_initial_state(population_states, populations, indices)
        @test length(initial) == 6
        result = first_order_correlation(model, populations, range(0.0, 1.0, length = 11))
        @test length(result.correlation) == 11

        stochastic = simulate_first_order_coherence_trajectory(
            model,
            2,
            2,
            range(0.0, 1.0, length = 11);
            rng = Xoshiro(1234),
        )
        @test length(stochastic) == 11
        @test stochastic[1] == A_JM_minus(2, 2)^2
    end

    @testset "full density representation" begin
        model = DickeModel(4; collective_decay = 1.0, local_decay = 0.1)
        space, jumps, rates, labels = collapse_operators(model)
        @test length(space) == 9
        @test length(jumps) == 4
        @test length(rates) == length(labels) == length(jumps)
        @test :collective_decay in labels
        @test :local_decay_minus in labels

        full = full_liouvillian(model)
        @test size(full.liouvillian.data) == (81, 81)
        initial = fully_excited_state(model)
        rho = dm(initial)
        derivative = reshape(full.liouvillian.data * vec(rho.data), size(rho.data))
        @test abs(tr(derivative)) < 1e-10
    end

    @testset "cavity composition" begin
        model = DickeModel(4; local_pump = 0.1)
        system = cavity_system(model; cutoff = 2, g = 0.1, kappa = 1.0)
        @test length(system.spin) == 9
        @test length(system.cavity) == 3
        @test length(system.basis) == 27
        @test last(system.labels) == :cavity_decay
        @test last(system.rates) == 1.0
    end

    @testset "spectra" begin
        times = range(0.0, 20.0, length = 2001)
        correlation = exp.(-0.5 .* times)
        frequencies, spectrum = spectrum_even_fft(times, correlation)
        @test length(frequencies) == length(spectrum)
        @test isapprox(spectrum_fwhm(frequencies, spectrum), 1.0; rtol = 0.05)
    end
end
