using LinearAlgebra
using Printf

function read_data_array(text, name)
    name_pos = findfirst("\"$name\"", text)
    name_pos === nothing && error("Could not find key '$name'")

    data_pos = findnext("\"data\"", text, last(name_pos))
    data_pos === nothing && error("Could not find data array for '$name'")

    open_pos = findnext('[', text, last(data_pos))
    close_pos = findnext(']', text, open_pos)
    open_pos === nothing && error("Could not find opening bracket for '$name'")
    close_pos === nothing && error("Could not find closing bracket for '$name'")

    entries = split(text[(open_pos + 1):(close_pos - 1)], ',')
    return [parse(Float64, strip(entry)) for entry in entries if !isempty(strip(entry))]
end

function read_scalar(text, name)
    name_pos = findfirst("\"$name\"", text)
    name_pos === nothing && error("Could not find key '$name'")

    data_pos = findnext("\"data\"", text, last(name_pos))
    data_pos === nothing && error("Could not find scalar data for '$name'")

    data_end = last(data_pos)
    value_match = match(r":\s*([-+0-9.eE]+)", text[data_end:end])
    value_match === nothing && error("Could not parse scalar '$name'")
    return parse(Int, value_match.captures[1])
end

function dynamics_matrix(n)
    A = zeros(2n, n)

    for t in 1:(2n), j in 1:n
        interval_start = 2j - 2
        interval_end = 2j

        if t <= interval_start
            A[t, j] = 0.0
        elseif t < interval_end
            A[t, j] = 0.5 * (t - interval_start)^2
        else
            A[t, j] = 2 * (t - 2j + 1)
        end
    end

    return A
end

function write_svg_plot(filename, coin_x, robot_x, missed_coin)
    width = 760
    height = 560
    margin = 55

    y = collect(1:length(coin_x))
    all_x = vcat(coin_x, robot_x)
    xmin, xmax = extrema(all_x)
    ymin, ymax = extrema(y)
    xpad = 0.08 * (xmax - xmin)
    ypad = 0.08 * (ymax - ymin)
    xmin -= xpad
    xmax += xpad
    ymin -= ypad
    ymax += ypad

    sx(x) = margin + (x - xmin) / (xmax - xmin) * (width - 2margin)
    sy(yval) = height - margin - (yval - ymin) / (ymax - ymin) * (height - 2margin)

    open(filename, "w") do io
        println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$width" height="$height" viewBox="0 0 $width $height">""")
        println(io, """<rect width="100%" height="100%" fill="white"/>""")
        println(io, """<line x1="$margin" y1="$(height - margin)" x2="$(width - margin)" y2="$(height - margin)" stroke="#222" stroke-width="1"/>""")
        println(io, """<line x1="$margin" y1="$margin" x2="$margin" y2="$(height - margin)" stroke="#222" stroke-width="1"/>""")
        println(io, """<text x="$(width / 2)" y="$(height - 15)" text-anchor="middle" font-family="sans-serif" font-size="14">x</text>""")
        println(io, """<text x="18" y="$(height / 2)" text-anchor="middle" font-family="sans-serif" font-size="14" transform="rotate(-90 18 $(height / 2))">y / time</text>""")

        print(io, """<polyline points=\"""")
        for (xval, yval) in zip(robot_x, y)
            @printf(io, "%.6f,%.6f ", sx(xval), sy(yval))
        end
        println(io, """" fill="none" stroke="#1f77b4" stroke-width="2.5"/>""")

        for (i, (xval, yval)) in enumerate(zip(coin_x, y))
            color = i == missed_coin ? "#d62728" : "#2ca02c"
            radius = i == missed_coin ? 5.0 : 4.0
            @printf(io, """<circle cx="%.6f" cy="%.6f" r="%.1f" fill="%s" opacity="0.9"/>\n""", sx(xval), sy(yval), radius, color)
        end

        for (xval, yval) in zip(robot_x, y)
            @printf(io, """<circle cx="%.6f" cy="%.6f" r="2.4" fill="#1f77b4"/>\n""", sx(xval), sy(yval))
        end

        println(io, """<text x="$(width - margin)" y="$(margin - 22)" text-anchor="end" font-family="sans-serif" font-size="14" fill="#1f77b4">robot path</text>""")
        println(io, """<text x="$(width - margin)" y="$(margin - 4)" text-anchor="end" font-family="sans-serif" font-size="14" fill="#2ca02c">collected coins</text>""")
        println(io, """<text x="$(width - margin)" y="$(margin + 14)" text-anchor="end" font-family="sans-serif" font-size="14" fill="#d62728">missed coin</text>""")
        println(io, "</svg>")
    end
end

function main()
    data_file = joinpath(@__DIR__, "robot_coin_collector.json")
    text = read(data_file, String)
    n = read_scalar(text, "n")
    coin_x = read_data_array(text, "x")

    A = dynamics_matrix(n)

    # 2(c)
    # println("Can collect all coins: ", is_consistent(A, coin_x))

    missed_coin = nothing
    best_f = nothing
    best_residual = nothing
    missed_error = nothing
    tol = 1e-8

    for miss in 1:(2n)
        rows = [i for i in 1:(2n) if i != miss]
        f = A[rows, :] \ coin_x[rows]

        residual = norm(A[rows, :] * f - coin_x[rows], Inf)

        if residual <= tol
            missed_coin = miss
            best_f = f
            best_residual = residual
            missed_error = abs((A * f - coin_x)[miss])
            break
        end
    end

    println("Missed coin: ", missed_coin)
    # println("Residual on collected coins: ", best_residual)
    println("Error at missed coin: ", missed_error)
    # 2(e)
    println("Input force: ", round.(best_f; digits=4))
    # println("Robot x positions: ", A * best_f)
    # println("Maximum error on collected coins: ",
    #         norm((A * best_f)[[i for i in 1:(2n) if i != missed_coin]] -
    #              coin_x[[i for i in 1:(2n) if i != missed_coin]], Inf))

    plot_file = joinpath(@__DIR__, "p2_robot_path.svg")
    write_svg_plot(plot_file, coin_x, A * best_f, missed_coin)
    # println("Plot written to ", plot_file)
end

main()
