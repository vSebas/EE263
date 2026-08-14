using LinearAlgebra
ENV["GKSwstype"] = "100"
using Plots
using Printf

function ellipse_points(A; npoints=500)
    E = eigen(Symmetric(A))
    theta = range(0, 2pi, length=npoints)
    circle = hcat(cos.(theta), sin.(theta))'
    points = E.vectors * Diagonal(1 ./ sqrt.(E.values)) * circle
    return points[1, :], points[2, :], E.values, E.vectors
end

function tick_locations(values; step=0.5)
    lower = step * floor(minimum(values) / step)
    upper = step * ceil(maximum(values) / step)
    return lower:step:upper
end

function plot_ellipse(A, filename; show_semiaxes=false, vectors=nothing,
                      fill_feasible=false, title="", tick_step=0.5)
    x, y, vals, vecs = ellipse_points(A)
    p = plot(x, y;
             seriestype=fill_feasible ? :shape : :path,
             aspect_ratio=:equal,
             xticks=tick_locations(x; step=tick_step),
             yticks=tick_locations(y; step=tick_step),
             xlabel="x1",
             ylabel="x2",
             fillalpha=fill_feasible ? 0.18 : 0.0,
             label=fill_feasible ? "feasible set" : "ellipse",
             linewidth=2,
             title=title)

    if show_semiaxes
        for i in 1:2
            axis = vecs[:, i] / sqrt(vals[i])
            plot!(p, [-axis[1], axis[1]], [-axis[2], axis[2]];
                  linestyle=:dash,
                  linewidth=2,
                  label=i == 1 ? "semiaxes" : false)
        end
    end

    if vectors !== nothing
        for i in 1:size(vectors, 2)
            b = vectors[:, i]
            plot!(p, [0.0, b[1]], [0.0, b[2]];
                  arrow=true,
                  linewidth=2,
                  label=i == 1 ? "sensor vectors" : false)
            annotate!(p, b[1], b[2], text("b$i", 8))
        end
    end

    savefig(p, filename)
    return p
end

function main()
    img_dir = normpath(joinpath(@__DIR__, "..", "latex", "img"))
    mkpath(img_dir)

    A_b = [1.0 0.0;
           0.0 2.0]
    plot_ellipse(A_b, joinpath(img_dir, "p4b_ellipse.pdf");
                 title="x' A x = 1")

    A_c = [0.2 -0.1;
           -0.1 0.4]
    _, _, vals_c, vecs_c = ellipse_points(A_c)
    plot_ellipse(A_c, joinpath(img_dir, "p4c_ellipse.pdf");
                 show_semiaxes=true,
                 title="x' A x = 1")

    b1 = [0.89, 0.45]
    b2 = [0.45, 0.89]
    b3 = [-0.71, 0.71]
    B = hcat(b1, b2, b3)
    A_d = B * B'
    _, _, vals_d, _ = ellipse_points(A_d)
    plot_ellipse(A_d, joinpath(img_dir, "p4d_sensors.pdf");
                 vectors=B,
                 fill_feasible=true,
                 title="norm(B' x) <= 1")

    println("part (c) eigenvalues = ")
    show(stdout, "text/plain", round.(vals_c; digits=6))
    println()
    println("part (c) semiaxis lengths = ")
    show(stdout, "text/plain", round.(1 ./ sqrt.(vals_c); digits=6))
    println()
    println("part (c) semiaxis directions = ")
    show(stdout, "text/plain", round.(vecs_c; digits=6))
    println()
    println("part (d) A = B*B' = ")
    show(stdout, "text/plain", round.(A_d; digits=6))
    println()
    println("part (d) semiaxis lengths = ")
    show(stdout, "text/plain", round.(1 ./ sqrt.(vals_d); digits=6))
    println()
end

main()
