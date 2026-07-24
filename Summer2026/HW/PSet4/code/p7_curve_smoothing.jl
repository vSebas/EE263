using LinearAlgebra
ENV["GKSwstype"] = "100"
using Plots
using Printf

include(joinpath(@__DIR__, "readclassjson.jl"))

function second_difference_matrix(n)
    D = zeros(n - 2, n)
    scale = n^2
    for i in 2:(n - 1)
        row = i - 1
        D[row, i - 1] = scale
        D[row, i] = -2scale
        D[row, i + 1] = scale
    end
    return D
end

function smooth_curve(f, D, mu)
    n = length(f)
    if isinf(mu)
        T = hcat(ones(n), collect(1:n) ./ n)
        return T * (T \ f)
    end

    H = Matrix{Float64}(I, n, n) / n + mu * (D' * D) / (n - 2)
    return H \ (f / n)
end

deviation(f, g) = sum((f .- g).^2) / length(f)
curvature(D, g) = sum((D * g).^2) / size(D, 1)

function main()
    data = readclassjson(joinpath(@__DIR__, "curve_smoothing.json"))
    n = data["n"]
    f = data["f"]
    D = second_difference_matrix(n)

    mu_values = 10.0 .^ range(-13, -4, length=120)
    g_values = [smooth_curve(f, D, mu) for mu in mu_values]
    d_values = [deviation(f, g) for g in g_values]
    c_values = [curvature(D, g) for g in g_values]

    g0 = smooth_curve(f, D, 0.0)
    ginf = smooth_curve(f, D, Inf)
    selected_mu = [1e-7, 1e-6, 1e-4]
    selected_g = [smooth_curve(f, D, mu) for mu in selected_mu]

    img_dir = normpath(joinpath(@__DIR__, "..", "latex", "img"))
    d_plot = vcat(deviation(f, g0), d_values, deviation(f, ginf))
    c_plot = vcat(curvature(D, g0), c_values, curvature(D, ginf))

    plt = plot(d_plot, c_plot;
               label="optimal tradeoff",
               xlabel="deviation d",
               ylabel="curvature c",
               linewidth=2)
    scatter!(plt, [d_plot[1]], [c_plot[1]];
             label="mu = 0", markersize=5)
    scatter!(plt, [d_plot[end]], [c_plot[end]];
             label="mu -> infinity", markersize=5)
    savefig(plt, joinpath(img_dir, "p7_tradeoff.pdf"))

    xs = collect(1:n)
    curve_plot = plot(xs, g0;
                      label="g, mu=0",
                      xlabel="index i",
                      ylabel="value",
                      linewidth=2,
                      legend=:topleft)
    for (mu, g) in zip(selected_mu, selected_g)
        plot!(curve_plot, xs, g; label="g, mu=$(mu)", linewidth=2)
    end
    plot!(curve_plot, xs, ginf; label="g, mu=infinity", linewidth=2)
    plot!(curve_plot, xs, f;
          label="data f",
          color=:black,
          linestyle=:dot,
          linewidth=2)
    savefig(curve_plot, joinpath(img_dir, "p7_smoothed_curves.pdf"))

    @printf("n = %d\n", n)
    @printf("mu = 0: d = %.8g, c = %.8g\n", deviation(f, g0), curvature(D, g0))
    for (mu, g) in zip(selected_mu, selected_g)
        @printf("mu = %.1e: d = %.8g, c = %.8g\n", mu, deviation(f, g), curvature(D, g))
    end
    @printf("mu -> infinity: d = %.8g, c = %.8g\n", deviation(f, ginf), curvature(D, ginf))
end

main()
