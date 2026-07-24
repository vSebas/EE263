function svg_to_pdf(svg_file, pdf_file)
    run(`rsvg-convert -f pdf -o $pdf_file $svg_file`)
end

function nice_bounds(values; pad_fraction=0.08)
    lo, hi = extrema(values)
    if lo == hi
        pad = max(abs(lo), 1.0) * pad_fraction
    else
        pad = (hi - lo) * pad_fraction
    end
    return lo - pad, hi + pad
end

function write_line_plot(filename, xs, series; labels=nothing, colors=nothing,
                         title="", xlabel="", ylabel="", width=780, height=520,
                         ylog=false, xlog=false, legend_x=nothing, legend_y=nothing)
    margin_left = 70
    margin_right = 25
    margin_top = 45
    margin_bottom = 60

    colors === nothing && (colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd",
                                    "#ff7f0e", "#17becf", "#8c564b"])
    labels === nothing && (labels = ["series $i" for i in 1:length(series)])

    tx(v) = xlog ? log10(v) : v
    ty(v) = ylog ? log10(max(v, eps(Float64))) : v

    all_x = [tx(x) for x in xs]
    all_y = Float64[]
    for ys in series
        append!(all_y, [ty(y) for y in ys])
    end

    xmin, xmax = nice_bounds(all_x)
    ymin, ymax = nice_bounds(all_y)

    sx(x) = margin_left + (tx(x) - xmin) / (xmax - xmin) * (width - margin_left - margin_right)
    sy(y) = height - margin_bottom - (ty(y) - ymin) / (ymax - ymin) * (height - margin_top - margin_bottom)

    open(filename, "w") do io
        println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$width" height="$height" viewBox="0 0 $width $height">""")
        println(io, """<rect width="100%" height="100%" fill="white"/>""")
        println(io, """<line x1="$margin_left" y1="$(height - margin_bottom)" x2="$(width - margin_right)" y2="$(height - margin_bottom)" stroke="#222" stroke-width="1"/>""")
        println(io, """<line x1="$margin_left" y1="$margin_top" x2="$margin_left" y2="$(height - margin_bottom)" stroke="#222" stroke-width="1"/>""")
        println(io, """<text x="$(width / 2)" y="25" text-anchor="middle" font-family="sans-serif" font-size="16">$title</text>""")
        println(io, """<text x="$(width / 2)" y="$(height - 18)" text-anchor="middle" font-family="sans-serif" font-size="13">$xlabel</text>""")
        println(io, """<text x="18" y="$(height / 2)" text-anchor="middle" font-family="sans-serif" font-size="13" transform="rotate(-90 18 $(height / 2))">$ylabel</text>""")

        for (sidx, ys) in enumerate(series)
            color = colors[mod1(sidx, length(colors))]
            dash = occursin("data", lowercase(labels[sidx])) ? " stroke-dasharray=\"4 4\"" : ""
            print(io, """<polyline points=\"""")
            for (x, y) in zip(xs, ys)
                print(io, "$(sx(x)),$(sy(y)) ")
            end
            println(io, """" fill="none" stroke="$color" stroke-width="2"$dash/>""")
        end

        legend_x === nothing && (legend_x = width - 190)
        legend_y === nothing && (legend_y = margin_top)
        for (sidx, label) in enumerate(labels)
            color = colors[mod1(sidx, length(colors))]
            y0 = legend_y + 18 * (sidx - 1)
            dash = occursin("data", lowercase(label)) ? " stroke-dasharray=\"4 4\"" : ""
            println(io, """<line x1="$legend_x" y1="$y0" x2="$(legend_x + 25)" y2="$y0" stroke="$color" stroke-width="2"$dash/>""")
            println(io, """<text x="$(legend_x + 32)" y="$(y0 + 4)" font-family="sans-serif" font-size="12">$label</text>""")
        end

        println(io, "</svg>")
    end
end

function write_heatmap_svg(filename, X; width=620, height=620)
    nrows, ncols = size(X)
    margin = 20
    xmin, xmax = extrema(vec(X))
    span = xmax - xmin
    span == 0 && (span = 1.0)
    cell_w = (width - 2margin) / ncols
    cell_h = (height - 2margin) / nrows

    gray(v) = begin
        q = clamp((v - xmin) / span, 0.0, 1.0)
        g = round(Int, 255 * (1 - q))
        "rgb($g,$g,$g)"
    end

    open(filename, "w") do io
        println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$width" height="$height" viewBox="0 0 $width $height">""")
        println(io, """<rect width="100%" height="100%" fill="white"/>""")
        for row in 1:nrows, col in 1:ncols
            x = margin + (col - 1) * cell_w
            y = margin + (row - 1) * cell_h
            println(io, """<rect x="$x" y="$y" width="$(cell_w + 0.2)" height="$(cell_h + 0.2)" fill="$(gray(X[row, col]))"/>""")
        end
        println(io, """<rect x="$margin" y="$margin" width="$(width - 2margin)" height="$(height - 2margin)" fill="none" stroke="#222" stroke-width="1"/>""")
        println(io, "</svg>")
    end
end
