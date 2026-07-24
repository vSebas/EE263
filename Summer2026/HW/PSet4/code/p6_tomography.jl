using LinearAlgebra
ENV["GKSwstype"] = "100"
using Plots
using Printf

include(joinpath(@__DIR__, "readclassjson.jl"))

function main()
    data = readclassjson(joinpath(@__DIR__, "tomo_data.json"))
    N = data["N"]
    npixels = data["npixels"]
    y = data["y"]
    line_pixel_lengths = data["line_pixel_lengths"]

    A = line_pixel_lengths'
    x = A \ y
    residual = A * x - y

    image = reshape(x, npixels, npixels)
    img_dir = normpath(joinpath(@__DIR__, "..", "latex", "img"))
    plt = heatmap(image;
                  yflip=true,
                  aspect_ratio=:equal,
                  color=:gist_gray,
                  cbar=:none,
                  framestyle=:none)
    savefig(plt, joinpath(img_dir, "p6_reconstruction.pdf"))

    @printf("npixels = %d\n", npixels)
    @printf("measurements = %d\n", N)
    @printf("unknowns = %d\n", npixels^2)
    @printf("least-squares residual norm = %.6g\n", norm(residual))
    @printf("RMS residual = %.6g\n", norm(residual) / sqrt(length(y)))
    @printf("density min = %.6g\n", minimum(x))
    @printf("density max = %.6g\n", maximum(x))
end

main()
