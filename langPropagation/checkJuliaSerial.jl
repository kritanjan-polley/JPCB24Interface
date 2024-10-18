start_time = time()
using Random, StatsBase, LinearAlgebra
using DelimitedFiles, Printf
using Distributed, Polynomials
using BSplineKit, NaNStatistics, QuadGK
include("functions.jl");


const mass :: Float64 = 3.0*16.0*1e-3;
const dvalue_liquid :: Float64 = 1.85e-1; 
const dvalue_gas :: Float64 = 1.4e3; ##
const factor :: Float64 = 4.814e-1;

const rt :: Float64 = 300.0*8.314; ## J/mol
const kBT :: Float64 = 0.5962; ## kcal/mol
const gamma_liquid :: Float64 = kBT/mass/dvalue_liquid*factor; ## 1/ps
const gamma_gas :: Float64 = kBT/mass/dvalue_gas*factor; ## 1/ps 
const kBTm :: Float64 = rt/mass*1e-4; ## A^2/ps^2
const sigma :: Float64 = sqrt(kBTm);## A/ps

potential_data :: Matrix{Float64} = readdlm("freefile.txt", comments=true);
friction_optimization :: Matrix{Float64} = readdlm("optimization_new.txt", comments=true);
friction_optimization = friction_optimization[sortperm(friction_optimization[:,1]),:];

potential_data_x = potential_data[:,1];
potential_data_y = potential_data[:,2]*factor;

width_x, shift_x, power = friction_optimization[2,2:end];
const high_x = gamma_liquid;
const cutoff_left = 15.0; const cutoff_right = 27.5;
const boxlen = 55.882675

println("width_x is : ", width_x, " shift_x is : ", shift_x, " high_x is : ",high_x, " tanh power is : ", power)
println("Gamma liquid phase is : ", gamma_liquid, " (1/ps)")
println("gamma gas phase is : ", gamma_gas, " (1/ps)")
println("Box length is : ", boxlen, " A" )

function add_count(x::Float64, region::String)
    if region == "left"
        return ifelse( abs(x) <= cutoff_left, 1.0, 0.0  )
    elseif region == "middle"
        return ifelse( cutoff_right > abs(x) > cutoff_left, 1.0, 0.0 )
    elseif region == "right"
        return ifelse( abs(x) >= cutoff_right, 1.0, 0.0  )
    else
        println("wrong input for region or invalid values of xn")
        exit()
    end
end;

function gamma_function(x::Float64)
    result = @. (ifelse(x > 0, (tanh(width_x * (-x + shift_x)) / 2 + 0.5)^power * (high_x - gamma_gas) + gamma_gas, (tanh(width_x * (x + shift_x)) / 2 + 0.5)^power * (high_x - gamma_gas) + gamma_gas))
    return result
end;

function gamma_function_deriv(x::Float64)
    res = ifelse( x>0, -0.5*power*width_x*(high_x-gamma_gas)*(sech( width_x*(shift_x-x) )^2)*(0.5+0.5*tanh(width_x*(shift_x-x)))^(power-1.0), 0.5*power*width_x*(high_x-gamma_gas)*(sech( width_x*(x + shift_x) )^2)*(0.5+0.5*tanh(width_x*(x + shift_x)))^(power-1.0)  )
    return res
end;

function gamma_function_primitive()
    ## get integral values
    xArr = range(-boxlen, boxlen, length=10^3)
    vals = zeros(length(xArr))
    for i = 1:length(xArr)
        value, _ = quadgk( x->gamma_function(x), 0, xArr[i])
        vals[i] = value
    end
    ## fit cubic spline
    fit_gamma_primiteve = extrapolate(interpolate(xArr, vals, BSplineOrder(6), Natural()), Flat());

    return fit_gamma_primiteve
end;
gamma_fit_primitive = gamma_function_primitive();

zdata = range(minimum(potential_data_x), maximum(potential_data_x), round(Int, 1e3));
fit_pmf = extrapolate(interpolate(potential_data_x, potential_data_y, BSplineOrder(6), Natural()), Flat());
fit_force = diff(fit_pmf, Derivative(1));

pot_interpolate = fit_pmf
function force(x::Float64)
    return -fit_force(x)
end;

temp = hcat( pot_interpolate.(zdata), zdata )[sortperm(hcat( pot_interpolate.(zdata), zdata )[:,1]),:][:,2];
minloc = temp[temp .>0.0][1]
println("Location of minima is at ", minloc)

function modify_array_for_binning(xleft::Float64, xright::Float64, xlen::Integer)
    x_temp = range(xleft, xright, xlen)
    dx_temp = x_temp[2] - x_temp[1]
    x_temp_2 = range( x_temp[1]- dx_temp/2, x_temp[end] + dx_temp/2, length(x_temp)+1 )
    return x_temp_2
end;

const x_abs = boxlen; const v_abs = 10.0; const x_array_len = 201; const v_array_len = 201;

# mybin_x = modify_array_for_binning(-x_abs, x_abs, x_array_len);
# mybin_v = modify_array_for_binning(-v_abs*sigma, v_abs*sigma, v_array_len);

const target_t = 2000.0; const dt = 0.004; const dt2 = dt*dt; const xlen = floor(Int, target_t/dt);
const num_cycle = round(Int, 1e3)

xn_in = rand(Normal(minloc, 0.5), num_cycle )
vn_in = -abs.(rand(Normal(0.0, sigma), num_cycle ))

println("Starting with propagation loops : ", time() - start_time, " seconds");

# time_check_array_position = round.(Int, range(1, 30, step = 1)./dt)
# time_check_array_velocity = round.(Int, range(1, 10, step = 1)./dt)

println("Number of elements in xarray : ", xlen)
println("Number of cycles : ", num_cycle)
println("Time step : ", dt*1e3, " fs")


# for i in time_check_array_position
#     if isfile("probability_position_"*string(round(Int, i*dt))*".txt") == true
#     else
#         file_x = open("probability_position_"*string(round(Int, i*dt))*".txt","w")
#         close(file_x)
#     end
# end;
# for i in time_check_array_velocity
#     if isfile("probability_velocity_"*string(round(Int, i*dt))*".txt") == true
#     else
#         file_x = open("probability_velocity_"*string(round(Int, i*dt))*".txt","w")
#         close(file_x)
#     end
# end;

# values_to_be_saved_position = zeros(length(time_check_array_position), num_cycle)
# values_to_be_saved_velocity = zeros(length(time_check_array_velocity), num_cycle)

global t_a = zeros(xlen + 1); 
global t_b = zeros(xlen + 1); 
global t_c = zeros(xlen + 1);

println("Starting with propagation loops : ", time() - start_time, " seconds");

function compute_hist()
    @fastmath @inbounds for k = 1:num_cycle
        xn = zeros(xlen + 1)
        vn = zeros(xlen + 1)

        xn[1] = xn_in[k] ## rand( Normal(1.0, 1.0)) ##xn_in[k]
        vn[1] = vn_in[k] ## rand( Normal(1.0, 0.5*sig)) ##vn_in[k]

        for i = 1:xlen
            # compute coefficients
            xn1_bar = xn[i] + dt*vn[i] ## + dt2*force(xn[i])*0.5/mass
            alpha_r = (gamma_fit_primitive(xn1_bar) - gamma_fit_primitive(xn[i]))/(xn1_bar - xn[i])
            alpha_t =  gamma_function(xn[i]) + gamma_function_deriv(xn[i])*dt*0.5*vn[i]
            b = 1.0/(1.0 + 0.5*dt*alpha_r)
            a = b*(1.0 - 0.5*dt*alpha_r)
            n1 = randn() ## rand(Normal(0.0,1.0))
            c = b*sqrt(abs(2*alpha_t*dt*kBTm))*n1

            # update position, force and velocity
            xn[i+1] = xn[i] + b*dt*vn[i] + b*0.5*force(xn[i])*dt2/mass + 0.5*c*dt
            # xn[i+1] = xn[i+1] - 0.25*gamma_function_deriv(xn[i])/gamma_function(xn[i])*(kBTm)*dt2
            vn[i+1] = a*vn[i] + 0.5*dt*(a*force(xn[i]) + force(xn[i+1]))/mass + c

            if xn[i+1] > boxlen
                xn[i+1] = -boxlen + mod(xn[i+1],boxlen)
            elseif xn[i+1] < -boxlen
                xn[i+1] = boxlen - mod(abs(xn[i+1]),boxlen)
            end

            # if i in time_check_array_position
            #     temp_int = round(Int, i*dt)
            #     values_to_be_saved_position[temp_int, k] = xn[i]
            # end

            # if i in time_check_array_velocity
            #     temp_int = round(Int, i*dt)
            #     values_to_be_saved_velocity[temp_int, k] = vn[i]
            # end
        end

        global t_a += map( x->add_count(x, "left"), xn)
        global t_b += map( x->add_count(x, "middle"), xn)
        global t_c += map( x->add_count(x, "right"), xn)

        if mod(k,num_cycle/20) == 0.0
            println("Done with ", k, " th iteration after ", time()-start_time," seconds")
        end
    end
    return [t_a/num_cycle, t_b/num_cycle, t_c/num_cycle]
end

t_a, t_b, t_c = compute_hist()

# ## position
# for i in time_check_array_position
#     temp_int = round(Int, i*dt)
#     file_x = open("probability_position_"*string(temp_int)*".txt","a")
#     writedlm(file_x, values_to_be_saved_position[temp_int, :])
#     close(file_x)
# end;

# ## velocity
# for i in time_check_array_velocity
#     temp_int = round(Int, i*dt)
#     file_v = open("probability_velocity_"*string(temp_int)*".txt","a")
#     writedlm(file_v, values_to_be_saved_velocity[temp_int, :])
#     close(file_v)
# end;
random_integer = round(Int, time_ns())
const n_every = 10

# writedlm("array_x_"*string(random_integer)*".txt",  hcat(moving_average(mybin_x,2), ax), ' ')
# writedlm("array_v_"*string(random_integer)*".txt",  hcat(moving_average(mybin_v,2), av), ' ')
writedlm("t_a_"*string(random_integer)*".txt",  hcat(range(0.0,target_t, xlen+1)[begin:n_every:end], t_a[begin:n_every:end]), ' ')
writedlm("t_b_"*string(random_integer)*".txt",  hcat(range(0.0,target_t, xlen+1)[begin:n_every:end], t_b[begin:n_every:end]), ' ')
writedlm("t_c_"*string(random_integer)*".txt",  hcat(range(0.0,target_t, xlen+1)[begin:n_every:end], t_c[begin:n_every:end]), ' ')

println("Done after ", time() - start_time," seconds")
