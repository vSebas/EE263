import JSON

function class_type(name)
    types = Dict(
        "Int64" => Int64,
        "Int32" => Int32,
        "Int8" => Int8,
        "UInt64" => UInt64,
        "UInt32" => UInt32,
        "UInt8" => UInt8,
        "Float64" => Float64,
        "Float32" => Float32,
        "String" => String,
        "Bool" => Bool,
    )
    return types[name]
end

function convert_class_data(data, ::Type{T}) where {T}
    if data isa Vector
        if !isempty(data) && first(data) isa Vector
            rows = length(data)
            cols = length(data[1])
            M = Matrix{T}(undef, rows, cols)
            for i in 1:rows, j in 1:cols
                M[i, j] = data[i][j]
            end
            return M
        end
        return convert(Vector{T}, data)
    end
    return convert(T, data)
end

function readclassjson(filename)
    raw = JSON.parsefile(filename)
    out = Dict{String, Any}()
    for name in keys(raw)
        T = class_type(raw[name]["type"])
        out[name] = convert_class_data(raw[name]["data"], T)
    end
    return out
end
