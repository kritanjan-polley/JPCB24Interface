using SparseArrays
using LinearAlgebra
# using SuiteSparseGraphBLAS
eps = 1e-9
large_eps = 1e-5


function find_nearest(array::AbstractArray, value::Float64)
    return argmin(abs.(array .- value))
end;

function eye(n::Integer, k::Integer)
    return diagm(k=>ones(n-abs(k)))
end;

deriv_accuracy = 2;

## 1st derivative
function deriv1(n::Integer, accuracy = deriv_accuracy, bc_left = "none", bc_right = "none")
    if accuracy == 2
        output = -0.5.*eye(n,-1) + 0.5.*eye(n,+1)
    elseif accuracy == 4
        output = 1/12.0.*eye(n,-2) -2/3.0.*eye(n,-1) +2/3.0.*eye(n,1) -1/12.0.*eye(n,2)
    elseif accuracy == 6
        output = -1/60*eye(n,-3)+3/20*eye(n,-2) -3/4*eye(n,-1) +3/4*eye(n,1) -3/20*eye(n,2)+1/60*eye(n,3)
    elseif accuracy == 8
        output = 1/280*eye(n,-4)-4/105*eye(n,-3)+1/5*eye(n,-2)-4/5*eye(n,-1)+4/5*eye(n,1)-1/5*eye(n,2)+4/105*eye(n,3)-1/280*eye(n,4)   
    end


    if bc_left == "periodic"
        if accuracy == 2
            output[1,n] = -0.5
        elseif accuracy == 4
            output[1,n] = -2/3
            output[1,n-1] = 1/12
            output[2,n] = 1/12
        elseif accuracy == 6
            output[1,n] = -3/4
            output[1,n-1] = 3/20
            output[1,n-2] = -1/60
            output[2,n] = 3/20
            output[2,n-1] = -1/60
            output[3,n] = 3/20
        elseif accuracy == 8
            output[1,n] = -4/5
            output[1,n-1] = 1/5
            output[1,n-2] = -4/105
            output[1,n-3] = 1/280
            output[2,n] = 1/5
            output[2,n-1] = -4/105
            output[2,n-2] = 1/280
            output[3,n] = -4/105
            output[3,n-1] = 1/280
            output[4,n] = 1/280
        end
    elseif bc_left == "reflect"
        if accuracy == 2
            output[1,2] = 0.0
        elseif accuracy == 4
            output[1,2] = 0.0
            output[1,3] = 0.0
            output[2,1] = 1/12
        elseif accuracy == 6
            output[1,2] = 0.0
            output[1,3] = 0.0
            output[1,4] = 0.0
            output[2,2] = 3/20 
            output[2,3] = 3/4 - 1/60
            output[3,2] = -1/60 - 3/4
        elseif accuracy == 8
            output[1,2] = 0.0
            output[1,3] = 0.0
            output[1,4] = 0.0
            output[1,5] = 0.0
            output[2,1] = 1/5
            output[2,2] = -4/105
            output[2,3] =  1/280
            output[3,1] = -4/105
            output[3,2] = 1/280
            output[4,1] = 1/280
        end
    elseif bc_left == "absorbing"
    elseif bc_left == "none"
    else
        println("Unknown left boundary condition")
        println("Exiting the program")
        exit()
    end

    if bc_right == "periodic"
        if accuracy == 2
            output[n,1] = 0.5
        elseif accuracy == 4
            println("not done yet")
            exit()
        elseif accuracy == 6
            println("not done yet")
            exit()
        elseif accuracy == 8
            println("not done yet")
            exit()
        end
    elseif bc_right == "reflect"
        if accuracy == 2
            output[n,n-1] = 0.0
        elseif accuracy == 4
            println("not done yet")
            exit()
        elseif accuracy == 6
            output[n,n-1] = 0.0
            output[n,n-2] = 0.0
            output[n,n-3] = 0.0
            output[n-1,n] = 3/20
            output[n-1,n-1] = -1/60
            output[n-1,n] = -1/60
        elseif accuracy == 8
            println("not done yet")
            exit()
        end
    elseif bc_right == "absorbing"
    elseif bc_right == "none"
    else
        println("Unknown right boundary condition")
        println("Exiting the program")
        exit()
    end

    return sparse(output)
end;

## 2nd derivative
function deriv2(n::Integer, accuracy = deriv_accuracy, bc_left = "none", bc_right = "none")
    if accuracy == 2
        factor = 1
        output = eye(n,-1)-2.0*eye(n,0)+eye(n,1)
    elseif accuracy == 4
        factor = 4/3
        output = -1/12*eye(n,-2)+4/3*eye(n,-1)-5/2*eye(n,0)+4/3*eye(n,1)-1/12*eye(n,2)
    elseif accuracy == 6
        factor = 3/2
        output = 1/90*eye(n,-3)-3/20*eye(n,-2)+3/2*eye(n,-1)-49/18*eye(n,0)+3/2*eye(n,1)-3/20*eye(n,2)+1/90*eye(n,3)
    elseif accuracy == 8
        factor = 8/5
        output = -1/560*eye(n,-4)+8/315*eye(n,-3)-1/5*eye(n,-2)+8/5*eye(n,-1)-205/72*eye(n,0)+8/5*eye(n,1)-1/5*eye(n,2)+8/315*eye(n,3)-1/560*eye(n,4)
    end

    if bc_left == "periodic"
        output[1,n] = 1.0
    elseif bc_left == "reflect"
        if accuracy == 2
            output[1,2] = 2.0
        elseif accuracy == 4
            println("not done yet")
            exit()
        elseif accuracy == 6
            output[1,2] = 2*3/2
            output[1,3] = 2*(-3/20)
            output[1,4] = 2*(1/90)
            output[2,2] = -49/18 - 3/20
            output[2,3] = 3/2 + 1/90
            output[3,2] = 3/2 + 1/90
        elseif accuracy == 8
            println("not done yet")
            exit()
        end
    elseif bc_left == "absorbing"
    elseif bc_left == "none"
    else
        println("Unknown left boundary condition")
        println("Exiting the program")
        exit()
    end

    if bc_right == "periodic"
        output[n,1] = 1.0
    elseif bc_right == "reflect"
        if accuracy == 2
            output[n,n-1] = 2.0
        elseif accuracy == 4
            output[n,n-1] = 2*4/3
            output[n,n-2] = 2*(-1/12)
            output[n-1,n] = 2*(-1/12)
        elseif accuracy == 6
            output[n,n-1] = 2*3/2
            output[n,n-2] = 2*(-3/20)
            output[n,n-3] = 2*(1/90)
            output[n-1,n] = -2*(3/20)
            output[n-1,n-1] = 2*(1/90)
            output[n-1,n] = 2*(1/90)
        elseif accuracy == 8
            println("not done yet")
            exit()
        end
    elseif bc_right == "absorbing"
    elseif bc_right == "none"
    else
        println("Unknown right boundary condition")
        println("Exiting the program")
        exit()
    end

    return sparse(output)
end;

## normal Distributions
function gauss(x::AbstractArray, mu::Float64, sig::Float64)
    return 1.0/sig/sqrt(2.0*pi).*exp.(-0.5*((x.-mu)./sig).^2) 
end;

function one_sided_gaussian(x::AbstractArray, mu::Float64, sig::Float64, side = "L")
    if side == "L"
        flip_heaviside = ifelse.( x.>= mu, 0.0, 1.0 )
        return 2.0 .* gauss(x, mu, sig).*flip_heaviside
    elseif side == "R"
        heaviside = ifelse.(x.>=mu, 1.0, 0.0)
        return  2.0 .* gauss(x, mu, sig).*heaviside
    end;
end;

## moving average 
function moving_average(x::AbstractArray, n::Integer)
    cum_sum = cumsum( append!([0.0], x) )
    return (cum_sum[n+1:1:end] - cum_sum[1:1:end-n])/float(n)
end;

function mychop(x::AbstractArray, mytol=eps)
    return map!(y -> isapprox(y, 0, atol=mytol) ? 0 : y, Array(x), Array(x))
end;
