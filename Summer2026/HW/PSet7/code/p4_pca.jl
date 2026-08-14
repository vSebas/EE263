using LinearAlgebra
using Statistics
ENV["GKSwstype"] = "100"
using Plots
using Printf

include(joinpath(@__DIR__, "readclassjson.jl"))

function digit_plot(x; title="", color=:grays, clims=nothing)
    image = reshape(x, 8, 8)'
    return heatmap(
        image;
        color,
        clims,
        yflip=true,
        aspect_ratio=:equal,
        axis=false,
        ticks=false,
        colorbar=false,
        title,
        titlefontsize=9,
    )
end

function main()
    data = readclassjson(joinpath(@__DIR__, "digits.json"))
    X = data["X"]
    n, N = size(X)

    xbar = mean(X; dims=2)
    Xc = X .- xbar
    F = svd(Xc)
    variance_fraction = cumsum(F.S .^ 2) ./ sum(F.S .^ 2)
    k90 = findfirst(>=(0.90), variance_fraction)

    img_dir = normpath(joinpath(@__DIR__, "..", "latex", "img"))
    mkpath(img_dir)

    p_variance = plot(
        1:n,
        variance_fraction;
        xlabel="number of principal components k",
        ylabel="fraction of variance captured",
        label="rho_k",
        linewidth=2,
        ylim=(0, 1.02),
        legend=:bottomright,
    )
    hline!(p_variance, [0.90]; linestyle=:dash, label="90%")
    scatter!(p_variance, [k90], [variance_fraction[k90]]; label="k = $k90")
    savefig(p_variance, joinpath(img_dir, "p4_variance.pdf"))

    ks = [2, 6, 12, 24]
    examples = [1, 2, 3]
    panels = Any[]
    for index in examples
        push!(panels, digit_plot(X[:, index]; title="original"))
        for k in ks
            Q = F.U[:, 1:k]
            xhat = vec(xbar) + Q * (Q' * Xc[:, index])
            push!(panels, digit_plot(xhat; title="k = $k"))
        end
    end
    p_recon = plot(panels...; layout=(length(examples), length(ks) + 1),
                   size=(900, 560))
    savefig(p_recon, joinpath(img_dir, "p4_reconstructions.pdf"))

    component_panels = Any[digit_plot(vec(xbar); title="mean")]
    for j in 1:5
        push!(component_panels,
              digit_plot(F.U[:, j]; title="u$j", color=:balance))
    end
    p_components = plot(component_panels...; layout=(1, 6), size=(900, 180))
    savefig(p_components, joinpath(img_dir, "p4_components.pdf"))

    k = 10
    Q = F.U[:, 1:k]
    Xhat = xbar .+ Q * (Q' * Xc)
    empirical_error = norm(X - Xhat)^2 / N
    predicted_error = sum(F.S[(k + 1):end] .^ 2) / N
    lost_fraction = 1 - variance_fraction[k]

    println("n = $n, N = $N")
    println("components for at least 90% variance = $k90")
    @printf("variance captured at k = 10: %.6f\n", variance_fraction[10])
    @printf("empirical mean-square error at k = 10: %.10f\n",
            empirical_error)
    @printf("discarded-singular-value prediction: %.10f\n",
            predicted_error)
    @printf("absolute verification difference: %.3e\n",
            abs(empirical_error - predicted_error))
    @printf("fraction of total variance lost: %.6f\n", lost_fraction)
end

main()
