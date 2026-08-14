using LinearAlgebra
ENV["GKSwstype"] = "100"
using Plots
using Printf

function terminal_matrices(T)
    Bp = [T - j + 0.5 for j in 1:T]
    Bv = ones(T)
    return Bp, Bv
end

function position_row(k, T)
    row = zeros(T)
    for j in 1:k
        row[j] = k - j + 0.5
    end
    return row
end

function trajectory(x; p0=0.0, v0=0.0, dt=0.02)
    times = collect(0.0:dt:length(x))
    p = similar(times)
    v = similar(times)
    f = similar(times)
    for (idx, t) in enumerate(times)
        if t == 0
            k = 0
            tau = 0.0
        else
            k = min(floor(Int, t), length(x))
            tau = t - k
            if tau == 0 && t > 0
                k -= 1
                tau = 1.0
            end
        end
        pcur = p0
        vcur = v0
        for j in 1:k
            pcur += vcur + 0.5 * x[j]
            vcur += x[j]
        end
        if k < length(x)
            pcur += vcur * tau + 0.5 * x[k + 1] * tau^2
            vcur += x[k + 1] * tau
            f[idx] = x[k + 1]
        else
            f[idx] = x[end]
        end
        p[idx] = pcur
        v[idx] = vcur
    end
    return times, p, v, f
end

function solve_part_a()
    T = 10
    p10, v10 = terminal_matrices(T)
    C = [p10'; v10'; position_row(5, T)']
    d = [1.0, 0.0, 0.0]
    x = C' * inv(C * C') * d
    return x, C * x
end

function solve_part_b(mu_values)
    T = 10
    p10, v10 = terminal_matrices(T)
    B = [p10'; v10']
    b = [10.0, 1.0]
    xs = [-(B' * B + mu * I) \ (B' * b) for mu in mu_values]
    J1 = [norm(B * x + b)^2 for x in xs]
    J2 = [sum(x.^2) for x in xs]
    x0 = -B' * inv(B * B') * b
    pushfirst!(xs, x0)
    pushfirst!(J1, norm(B * x0 + b)^2)
    pushfirst!(J2, sum(x0.^2))
    push!(xs, zeros(T))
    push!(J1, norm(b)^2)
    push!(J2, 0.0)
    return xs, J1, J2
end

function main()
    img_dir = normpath(joinpath(@__DIR__, "..", "latex", "img"))

    x_a, checks = solve_part_a()
    times, p, v, f = trajectory(x_a)
    force_plot = plot(times, f;
                      seriestype=:steppost,
                      label="f(t)",
                      xlabel="t",
                      ylabel="force")
    savefig(force_plot, joinpath(img_dir, "p4a_force.pdf"))

    state_plot = plot(times, p;
                      label="p(t)",
                      xlabel="t",
                      ylabel="state")
    plot!(state_plot, times, v; label="pdot(t)")
    savefig(state_plot, joinpath(img_dir, "p4a_state.pdf"))

    mu_values = 10.0 .^ range(-6, 6, length=50)
    _, J1, J2 = solve_part_b(mu_values)
    trade_plot = plot(J2, J1;
                      label="optimal tradeoff",
                      xlabel="J2",
                      ylabel="J1",
                      linewidth=2)
    scatter!(trade_plot, [J2[1]], [J1[1]];
             label="J1 = 0 endpoint",
             markersize=5)
    scatter!(trade_plot, [J2[end]], [J1[end]];
             label="J2 = 0 endpoint",
             markersize=5)
    savefig(trade_plot, joinpath(img_dir, "p4b_tradeoff.pdf"))

    println("part (a) x = ")
    show(stdout, "text/plain", round.(x_a; digits=6))
    println()
    println("[p(10), pdot(10), p(5)] = ")
    show(stdout, "text/plain", round.(checks; digits=6))
    println()
    @printf("part (a) energy = %.6f\n", sum(x_a.^2))
    @printf("part (b) J1=0 endpoint: J1 = %.6g, J2 = %.6f\n", J1[1], J2[1])
    @printf("part (b) J2=0 endpoint: J1 = %.6f, J2 = %.6g\n", J1[end], J2[end])
end

main()
