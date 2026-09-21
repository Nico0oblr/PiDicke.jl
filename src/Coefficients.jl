# These amplitudes are taken directly from the existing PermRates.jl and
# FullPerm.jl implementations.

A_JM_minus(J, M) = sqrt((J + M) * (J - M + 1))
A_JM_plus(J, M) = sqrt((J - M) * (J + M + 1))

function P_JM_minus_0(J, M, N)
    J == 0 && return 0.0
    return sqrt((N + 2) / (4 * J * (J + 1))) * A_JM_minus(J, M)
end

function P_JM_minus_minus(J, M, N)
    J == 0 && return 0.0
    return -sqrt((N + 2 * J + 2) * (J + M) * (J + M - 1) / (4 * J * (2 * J + 1)))
end

P_JM_minus_plus(J, M, N) =
    sqrt((N - 2 * J) * (J - M + 1) * (J - M + 2) / (4 * (J + 1) * (2 * J + 1)))

function P_JM_plus_0(J, M, N)
    J == 0 && return 0.0
    return sqrt((N + 2) / (4 * J * (J + 1))) * A_JM_plus(J, M)
end

function P_JM_plus_minus(J, M, N)
    J == 0 && return 0.0
    return sqrt((N + 2 * J + 2) * (J - M) * (J - M - 1) / (4 * J * (2 * J + 1)))
end

P_JM_plus_plus(J, M, N) =
    -sqrt((N - 2 * J) * (J + M + 1) * (J + M + 2) / (4 * (J + 1) * (2 * J + 1)))

function P_JM_z_0(J, M, N)
    J == 0 && return 0.0
    return sqrt((N + 2) / (4 * J * (J + 1))) * M
end

function P_JM_z_minus(J, M, N)
    J == 0 && return 0.0
    return sqrt((N + 2 * J + 2) * (J - M) * (J + M) / (4 * J * (2 * J + 1)))
end

P_JM_z_plus(J, M, N) =
    sqrt((N - 2 * J) * (J + 1 - M) * (J + 1 + M) / (4 * (J + 1) * (2 * J + 1)))

A_JM_minus2(J, M) = (J + M) * (J - M + 1)
A_JM_plus2(J, M) = (J - M) * (J + M + 1)

function P_JM_minus_02(J, M, N)
    J == 0 && return 0.0
    return (N + 2) * A_JM_minus2(J, M) / (4 * J * (J + 1))
end

function P_JM_minus_minus2(J, M, N)
    J == 0 && return 0.0
    return (N + 2 * J + 2) * (J + M) * (J + M - 1) / (4 * J * (2 * J + 1))
end

P_JM_minus_plus2(J, M, N) =
    (N - 2 * J) * (J - M + 1) * (J - M + 2) / (4 * (J + 1) * (2 * J + 1))

function P_JM_plus_02(J, M, N)
    J == 0 && return 0.0
    return (N + 2) * A_JM_plus2(J, M) / (4 * J * (J + 1))
end

function P_JM_plus_minus2(J, M, N)
    J == 0 && return 0.0
    return (N + 2 * J + 2) * (J - M) * (J - M - 1) / (4 * J * (2 * J + 1))
end

P_JM_plus_plus2(J, M, N) =
    (N - 2 * J) * (J + M + 1) * (J + M + 2) / (4 * (J + 1) * (2 * J + 1))

function P_JM_z_02(J, M, N)
    J == 0 && return 0.0
    return (N + 2) * M^2 / (4 * J * (J + 1))
end

function P_JM_z_minus2(J, M, N)
    J == 0 && return 0.0
    return (N + 2 * J + 2) * (J - M) * (J + M) / (4 * J * (2 * J + 1))
end

P_JM_z_plus2(J, M, N) =
    (N - 2 * J) * (J + 1 - M) * (J + 1 + M) / (4 * (J + 1) * (2 * J + 1))
