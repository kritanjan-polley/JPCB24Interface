using SparseArrays
using LinearAlgebra
using Random, Distributions
const eps = 1e-9

find_nearest(array :: AbstractArray, values :: Float64) = argmin(abs.(array .- value))
eye(n::Int64, k :: Int64) = diagm(k=>ones(n-abs(k)))

deriv_accuracy = 2; x_periodicity = true; v_periodicity = true;

## 1st derivative
function deriv1(n, accuracy = deriv_accuracy, periodic = true)
    if accuracy == 2
        output = -0.5.*eye(n,-1) + 0.5.*eye(n,+1)
    elseif accuracy == 4
        output = 1/12.0.*eye(n,-2) -2/3.0.*eye(n,-1) +2/3.0.*eye(n,1) -1/12.0.*eye(n,2)
    elseif accuracy == 6
        output = -1/60*eye(n,-3)+3/20*eye(n,-2) -3/4*eye(n,-1) +3/4*eye(n,1) -3/20*eye(n,2)+1/60*eye(n,3)
    elseif accuracy == 8
        output = 1/280*eye(n,-4)-4/105*eye(n,-3)+1/5*eye(n,-2)-4/5*eye(n,-1)+4/5*eye(n,1)-1/5*eye(n,2)+4/105*eye(n,3)-1/280*eye(n,4)   
    end

    if periodic == true
        output[1,n] = 1.0
        output[n,1] = -1.0
    end
    return sparse(output)
end;

## 2nd derivative
function deriv2(n, accuracy = deriv_accuracy, periodic = true)
    if accuracy == 2
        output = eye(n,-1)-2.0*eye(n,0)+eye(n,1)
    elseif accuracy == 4
        output = -1/12*eye(n,-2)+4/3*eye(n,-1)-5/2*eye(n,0)+4/3*eye(n,1)-1/12*eye(n,2)
    elseif accuracy == 6
        output = 1/90*eye(n,-3)-3/20*eye(n,-2)+3/2*eye(n,-1)-49/18*eye(n,0)+3/2*eye(n,1)-3/20*eye(n,2)+1/90*eye(n,3)
    elseif accuracy == 8
        output = -1/560*eye(n,-4)+8/315*eye(n,-3)-1/5*eye(n,-2)+8/5*eye(n,-1)-205/72*eye(n,0)+8/5*eye(n,1)-1/5*eye(n,2)+8/315*eye(n,3)-1/560*eye(n,4)
    end

    if periodic == true
        output[1,n] = 1.0
        output[n,1] = 1.0
    end

    return sparse(output)
end;

## normal Distribution
function gauss(x, mu, sig)
    return 1.0/sig/sqrt(2.0*pi).*exp.(-0.5*((x.-mu)./sig).^2) 
end;

## mathematica chop function
function mychop(x, mytol=eps)
    return map!(y -> isapprox(y, 0, atol=mytol) ? 0 : y, collect(x), collect(x))
end;

## moving average
function moving_average(x, n)
    cum_sum = cumsum( append!([0.0], x) )
    return (cum_sum[n+1:1:end] - cum_sum[1:1:end-n])/float(n)
end;

function wrapArray(x, left_val, right_val)
    return_x = similar(x)
    boxLen = abs(left_val) + abs(right_val)
    return_x = ifelse.( x.<= left_val, boxLen .- abs.(x), ifelse.( x.>=right_val, -boxLen .+x, x) )
    return return_x
end;

# function myDeriv(ys, xs)
#     fi = linear_interpolation(xs, ys, extrapolation_bc=Flat())
#     fid(x) = only(Interpolations.gradient(fi, x))
#     return fid.(xs)
# end;