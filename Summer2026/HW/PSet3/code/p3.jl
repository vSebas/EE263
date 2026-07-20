using LinearAlgebra
using Printf

function solve_medium(t, r; x1=1.0)
    n = length(t) + 1
    nvars = 2n - 1
    A = zeros(nvars, nvars)
    b = zeros(nvars)

    idx_x(i) = i - 1
    idx_y(i) = (n - 1) + i

    row = 1

    for i in 1:n-1
        # x_{i+1} - t_ix_i - r_iy_{i+1} = 0
        A[row, idx_x(i + 1)] = 1.0
        if i == 1
            b[row] = t[i] * x1
        else
            A[row, idx_x(i)] = -t[i]
        end
        A[row, idx_y(i + 1)] = -r[i]
        row += 1

        # y_{i} - r_ix_i - t_iy_{i+1} = 0
        A[row, idx_y(i)] = 1.0
        if i == 1
            b[row] = r[i] * x1
        else
            A[row, idx_x(i)] = -r[i]
        end
        A[row, idx_y(i + 1)] = -t[i]
        row += 1
    end

    # Boundary condition: y_{n} = x_{n}
    A[row, idx_y(n)] = 1.0
    A[row, idx_x(n)] = -1.0

    z = A \ b

    x = [x1; z[1:(n - 1)]]
    y = z[n:end]

    return x, y, y[1] / x1
end

function write_svg_plot(filename, x, y)
    width = 780
    height = 520
    margin = 55
    layers = collect(1:length(x))

    ymin, ymax = extrema(vcat(x, y))
    ypad = 0.08 * (ymax - ymin)
    ymin -= ypad
    ymax += ypad

    sx(i) = margin + (i - 1) / (length(x) - 1) * (width - 2margin)
    sy(v) = height - margin - (v - ymin) / (ymax - ymin) * (height - 2margin)

    open(filename, "w") do io
        println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$width" height="$height" viewBox="0 0 $width $height">""")
        println(io, """<rect width="100%" height="100%" fill="white"/>""")
        println(io, """<line x1="$margin" y1="$(height - margin)" x2="$(width - margin)" y2="$(height - margin)" stroke="#222" stroke-width="1"/>""")
        println(io, """<line x1="$margin" y1="$margin" x2="$margin" y2="$(height - margin)" stroke="#222" stroke-width="1"/>""")
        println(io, """<text x="$(width / 2)" y="$(height - 15)" text-anchor="middle" font-family="sans-serif" font-size="14">layer i</text>""")
        println(io, """<text x="18" y="$(height / 2)" text-anchor="middle" font-family="sans-serif" font-size="14" transform="rotate(-90 18 $(height / 2))">amplitude</text>""")

        print(io, """<polyline points=\"""")
        for (i, v) in zip(layers, x)
            @printf(io, "%.6f,%.6f ", sx(i), sy(v))
        end
        println(io, """" fill="none" stroke="#1f77b4" stroke-width="2.5"/>""")

        print(io, """<polyline points=\"""")
        for (i, v) in zip(layers, y)
            @printf(io, "%.6f,%.6f ", sx(i), sy(v))
        end
        println(io, """" fill="none" stroke="#d62728" stroke-width="2.5"/>""")

        for (i, v) in zip(layers, x)
            @printf(io, """<circle cx="%.6f" cy="%.6f" r="3" fill="#1f77b4"/>\n""", sx(i), sy(v))
        end
        for (i, v) in zip(layers, y)
            @printf(io, """<circle cx="%.6f" cy="%.6f" r="3" fill="#d62728"/>\n""", sx(i), sy(v))
        end

        println(io, """<text x="$(width - margin)" y="$(margin - 18)" text-anchor="end" font-family="sans-serif" font-size="14" fill="#1f77b4">right-traveling x_i</text>""")
        println(io, """<text x="$(width - margin)" y="$(margin + 2)" text-anchor="end" font-family="sans-serif" font-size="14" fill="#d62728">left-traveling y_i</text>""")
        println(io, "</svg>")
    end
end

function main()
    n = 20
    t = fill(0.96, n - 1)
    r = fill(0.02, n - 1)

    x, y, S = solve_medium(t, r)
    println("Nominal scattering coefficient S = ", S)

    target = 0.70
    best_k = 0
    best_S = 0.0
    best_err = Inf

    for k in 1:(n - 1)
        tk = copy(t)
        rk = copy(r)
        tk[k] = 0.02
        rk[k] = 0.96
        _, _, Sk = solve_medium(tk, rk)
        err = abs(Sk - target)

        if err < best_err
            best_k = k
            best_S = Sk
            best_err = err
        end
    end

    println("S_target = ", target)
    println("k = ", best_k)
    println("Predicted S  = ", best_S)
    println("Error = ", best_err)

    plot_file = joinpath(@__DIR__, "p3_layered_medium.svg")
    write_svg_plot(plot_file, x, y)
    println("Plot written to ", plot_file)
end

main()
