"""Cosine-transform a one-sided correlation function on an arbitrary frequency grid."""
function spectrum_from_correlator(times, correlation, frequencies; damping::Real = 0.0)
    length(times) >= 2 || throw(ArgumentError("at least two time points are required"))
    dt = times[2] - times[1]
    spectrum = zeros(Float64, length(frequencies))
    @inbounds for i in eachindex(frequencies)
        frequency = frequencies[i]
        value = 0.0
        for j in eachindex(times)
            value += real(correlation[j]) * exp(-damping * times[j]) * cos(frequency * times[j])
        end
        spectrum[i] = 2 * dt * value
    end
    return spectrum
end

"""FFT of the even extension of a correlation function sampled from zero."""
function spectrum_even_fft(times, correlation; damping::Real = 0.0)
    N = length(times)
    N >= 2 || throw(ArgumentError("at least two time points are required"))
    dt = times[2] - times[1]
    all(isapprox(times[j + 1] - times[j], dt; rtol = 1e-10, atol = 1e-12) for j in 1:(N - 1)) ||
        throw(ArgumentError("times must be uniformly spaced"))
    isapprox(first(times), 0; atol = 1e-12) || throw(ArgumentError("times must start at zero"))

    damped = correlation .* exp.(-damping .* times)
    extended = vcat(damped, damped[(end - 1):-1:2])
    transform = fft(extended) * dt
    count = length(extended)
    frequencies = 2π .* vcat(0:(count ÷ 2), -((count - 1) ÷ 2):-1) ./ (count * dt)
    permutation = sortperm(frequencies)
    return frequencies[permutation], real.(transform[permutation])
end

"""Full width at half maximum of a sampled, single-peaked spectrum."""
function spectrum_fwhm(frequencies, spectrum)
    normalized = spectrum ./ maximum(spectrum)
    peak = argmax(normalized)
    left = findlast(i -> normalized[i] <= 0.5 <= normalized[i + 1], 1:(peak - 1))
    right_offset = findfirst(i -> normalized[i] >= 0.5 >= normalized[i + 1], peak:(length(normalized) - 1))
    isnothing(left) && return NaN
    isnothing(right_offset) && return NaN
    right = peak + right_offset - 1
    omega_left = frequencies[left] +
        (0.5 - normalized[left]) * (frequencies[left + 1] - frequencies[left]) /
        (normalized[left + 1] - normalized[left])
    omega_right = frequencies[right] +
        (0.5 - normalized[right]) * (frequencies[right + 1] - frequencies[right]) /
        (normalized[right + 1] - normalized[right])
    return omega_right - omega_left
end
