using PIDicke

model = DickeModel(
    20;
    collective_decay = 1.0,
    local_pump = 0.4,
    local_dephasing = 0.02,
)

steady = population_steady_state(model)
times = range(0.0, 20.0, length = 2001)
result = first_order_correlation(model, steady, times)
frequencies, spectrum = spectrum_even_fft(times, result.correlation)

println("linewidth = ", spectrum_fwhm(frequencies, spectrum))
