using LinearAlgebra
using Printf
using Statistics

include(joinpath(@__DIR__, "readclassjson.jl"))
include(joinpath(@__DIR__, "svg_utils.jl"))

function main()
    data = readclassjson(joinpath(@__DIR__, "inductor_data.json"))
    nturns = data["n"]
    w = data["w"]
    d = data["d"]
    D = data["D"]
    L = data["L"]

    X = hcat(ones(length(L)), log.(nturns), log.(w), log.(d), log.(D))
    theta = X \ log.(L)
    alpha = exp(theta[1])
    beta = theta[2:end]

    Lhat = exp.(X * theta)
    pct_error = 100 .* abs.(Lhat .- L) ./ L

    img_dir = normpath(joinpath(@__DIR__, "..", "latex", "img"))
    write_line_plot(joinpath(img_dir, "p5_inductor_fit.svg"),
                    collect(1:length(L)),
                    [L, Lhat];
                    labels=["measured data", "model"],
                    title="Measured and fitted inductances",
                    xlabel="inductor index", ylabel="inductance (nH)")
    svg_to_pdf(joinpath(img_dir, "p5_inductor_fit.svg"),
               joinpath(img_dir, "p5_inductor_fit.pdf"))

    @printf("alpha = %.8g\n", alpha)
    for i in eachindex(beta)
        @printf("beta%d = %.8g\n", i, beta[i])
    end
    @printf("average percentage error = %.4f%%\n", mean(pct_error))
    @printf("maximum percentage error = %.4f%%\n", maximum(pct_error))
    @printf("RMS log residual = %.6g\n", norm(X * theta - log.(L)) / sqrt(length(L)))
end

main()
