using LinearAlgebra

A = [-1 0 0 -1 1;
      0 1 1 0 0;
      1 0 0 1 0]

I3 = Matrix{Rational{Int}}(I, 3, 3)

function is_lower_triangular_rectangular(B)
    for i in axes(B, 1), j in axes(B, 2)
        if i < j && B[i, j] != 0
            return false
        end
    end
    return true
end

function verify_case(label, B, property_name, property_holds)
    println("Case $label")
    println("B =")
    display(B)
    println("A * B =")
    display(A * B)
    println("AB == I: ", A * B == I3)
    println(property_name, ": ", property_holds)
    println()
end

# (a) The second row of B is zero.
B_a = [0 0 0;
       0 0 0;
       0 1 0;
       0 0 1;
       1 0 1]

verify_case(
    "(a)",
    B_a,
    "second row is zero",
    all(B_a[2, :] .== 0),
)

# (d) The second and third rows of B are the same.
B_d = [0 0 0;
       0 1//2 0;
       0 1//2 0;
       0 0 1;
       1 0 1]

verify_case(
    "(d)",
    B_d,
    "second and third rows are equal",
    B_d[2, :] == B_d[3, :],
)

# (f) B is lower triangular, i.e., B_ij = 0 for i < j.
B_f = B_a

verify_case(
    "(f)",
    B_f,
    "lower triangular",
    is_lower_triangular_rectangular(B_f),
)
