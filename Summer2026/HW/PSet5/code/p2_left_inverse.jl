using LinearAlgebra

function main()
    G = [2//1 3//1;
         1//1 0//1;
         0//1 4//1;
         1//1 1//1;
        -1//1 2//1]

    Gtilde = [-3//1 -1//1;
              -1//1  0//1;
               2//1 -3//1;
              -1//1 -3//1;
               1//1  2//1]

    H = [-18//25  71//25 16//25 -10//25 0//1;
           8//25  -1//25  4//25 -15//25 0//1]

    println("H = ")
    show(stdout, "text/plain", H)
    println()
    println("H*G = ")
    show(stdout, "text/plain", H * G)
    println()
    println("H*Gtilde = ")
    show(stdout, "text/plain", H * Gtilde)
    println()
end

main()
