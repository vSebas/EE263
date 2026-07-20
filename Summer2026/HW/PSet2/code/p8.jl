using LinearAlgebra
using JSON3
using Printf

function read_interpolation_data(filename)
    data = JSON3.read(read(filename, String))
    x = Float64.(data.x.data)
    y = Float64.(data.y.data)
    return x, y
end

function interpolation_matrix(x, y, m)
    numerator_terms = [x .^ j for j in 0:m]
    denominator_terms = [-(y .* (x .^ j)) for j in 1:m]
    return hcat(numerator_terms..., denominator_terms...)
end

function rational_value(x, a, b)
    numerator = sum(a[j + 1] * x^j for j in 0:(length(a) - 1))
    denominator = 1.0 + sum(b[j] * x^j for j in 1:length(b))
    return numerator / denominator
end

function fit_smallest_degree(x, y; tol=1e-8, max_degree=div(length(x) - 1, 2))
    for m in 0:max_degree
        A = interpolation_matrix(x, y, m)
        theta = A \ y
        residual = norm(A * theta - y, Inf)

        if residual <= tol
            a = theta[1:(m + 1)]
            b = theta[(m + 2):end]
            return m, a, b, residual
        end
    end

    error("No interpolating rational function found up to degree $max_degree")
end

function write_svg_plot(filename, x, y, a, b)
    width = 900
    height = 560
    margin = 55

    xmin, xmax = extrema(x)
    xs = collect(range(xmin, xmax, length=1200))
    ys = [rational_value(t, a, b) for t in xs]

    finite_ys = [v for v in vcat(y, ys) if isfinite(v)]
    ymin, ymax = extrema(finite_ys)
    ypad = 0.08 * (ymax - ymin)
    ymin -= ypad
    ymax += ypad

    sx(t) = margin + (t - xmin) / (xmax - xmin) * (width - 2margin)
    sy(t) = height - margin - (t - ymin) / (ymax - ymin) * (height - 2margin)

    open(filename, "w") do io
        println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$width" height="$height" viewBox="0 0 $width $height">""")
        println(io, """<rect width="100%" height="100%" fill="white"/>""")
        println(io, """<line x1="$margin" y1="$(height - margin)" x2="$(width - margin)" y2="$(height - margin)" stroke="#222" stroke-width="1"/>""")
        println(io, """<line x1="$margin" y1="$margin" x2="$margin" y2="$(height - margin)" stroke="#222" stroke-width="1"/>""")
        println(io, """<text x="$(width / 2)" y="$(height - 16)" text-anchor="middle" font-family="sans-serif" font-size="14">x</text>""")
        println(io, """<text x="18" y="$(height / 2)" text-anchor="middle" font-family="sans-serif" font-size="14" transform="rotate(-90 18 $(height / 2))">y</text>""")

        path_started = false
        print(io, """<path d=\"""")
        for (t, v) in zip(xs, ys)
            if isfinite(v)
                cmd = path_started ? "L" : "M"
                @printf(io, "%s %.6f %.6f ", cmd, sx(t), sy(v))
                path_started = true
            else
                path_started = false
            end
        end
        println(io, """" fill="none" stroke="#1f77b4" stroke-width="2.5"/>""")

        for (xi, yi) in zip(x, y)
            @printf(io, """<circle cx="%.6f" cy="%.6f" r="3.2" fill="#d62728" opacity="0.85"/>\n""", sx(xi), sy(yi))
        end

        println(io, """<text x="$(width - margin)" y="$(margin - 18)" text-anchor="end" font-family="sans-serif" font-size="14" fill="#1f77b4">rational fit</text>""")
        println(io, """<text x="$(width - margin)" y="$(margin + 2)" text-anchor="end" font-family="sans-serif" font-size="14" fill="#d62728">data</text>""")
        println(io, "</svg>")
    end
end

function print_coefficients(label, coeffs; first_index=0)
    println(label)
    for (i, value) in enumerate(coeffs)
        @printf("  %s_%d = %.15g\n", label, i + first_index - 1, value)
    end
end

function main()
    data_file = joinpath(@__DIR__, "rational_interpolation_data.json")
    x, y = read_interpolation_data(data_file)

    m, a, b, residual = fit_smallest_degree(x, y)
    yhat = [rational_value(xi, a, b) for xi in x]
    max_error = norm(yhat - y, Inf)

    println("Smallest degree m = ", m)
    print_coefficients("a", a)
    print_coefficients("b", b; first_index=1)
    # @printf("linear-system residual max norm = %.3e\n", residual)
    # @printf("interpolation max absolute error = %.3e\n", max_error)
    @printf("error = %.3e\n", max_error)

    plot_file = joinpath(@__DIR__, "p8_rational_interpolation.svg")
    write_svg_plot(plot_file, x, y, a, b)
    println("Plot written to ", plot_file)
end

main()
