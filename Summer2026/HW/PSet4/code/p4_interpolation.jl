using LinearAlgebra

include(joinpath(@__DIR__, "readclassjson.jl"))
include(joinpath(@__DIR__, "svg_utils.jl"))

function linear_interpolation_matrix(n, sample_times)
    m = length(sample_times)
    A = zeros(n, m)

    for segment in 1:(m - 1)
        t_left = sample_times[segment]
        t_right = sample_times[segment + 1]

        for t in t_left:t_right
            alpha = (t - t_left) / (t_right - t_left)
            A[t, segment] = 1 - alpha
            A[t, segment + 1] = alpha
        end
    end

    return A
end

function select_independent_rows(B; tol=1e-10)
    selected = Int[]
    current = zeros(0, size(B, 2))
    current_rank = 0

    for i in 1:size(B, 1)
        candidate = vcat(current, reshape(B[i, :], 1, :))
        candidate_rank = rank(candidate; atol=tol, rtol=0.0)

        if candidate_rank > current_rank
            push!(selected, i)
            current = candidate
            current_rank = candidate_rank
        end

        current_rank == size(B, 2) && break
    end

    return selected
end

function main()
    n = 20
    sample_times = [1, 5, 11, 20]
    A = linear_interpolation_matrix(n, sample_times)

    img_dir = normpath(joinpath(@__DIR__, "..", "latex", "img"))
    write_line_plot(joinpath(img_dir, "p4_linear_columns.svg"),
                    collect(1:n),
                    [A[:, j] for j in 1:size(A, 2)];
                    labels=["a1", "a2", "a3", "a4"],
                    title="Linear interpolation basis",
                    xlabel="index", ylabel="value")
    svg_to_pdf(joinpath(img_dir, "p4_linear_columns.svg"),
               joinpath(img_dir, "p4_linear_columns.pdf"))

    data = readclassjson(joinpath(@__DIR__, "interp.json"))
    B = data["B"]
    t_special = select_independent_rows(B)
    Z = B[t_special, :]
    V = B * inv(Z)

    write_line_plot(joinpath(img_dir, "p4_special_basis.svg"),
                    collect(1:size(B, 1)),
                    [V[:, j] for j in 1:size(V, 2)];
                    labels=["v1", "v2", "v3"],
                    title="Special basis from interp.json",
                    xlabel="index", ylabel="value")
    svg_to_pdf(joinpath(img_dir, "p4_special_basis.svg"),
               joinpath(img_dir, "p4_special_basis.pdf"))

    println("Linear interpolation A:")
    show(stdout, "text/plain", round.(A; digits=4))
    println()
    println("rank(A) = ", rank(A))
    println("Selected sample times for interp.json: ", t_special)
    println("Z = ")
    show(stdout, "text/plain", round.(Z; digits=4))
    println()
    println("V = ")
    show(stdout, "text/plain", round.(V; digits=4))
    println()
end

main()
