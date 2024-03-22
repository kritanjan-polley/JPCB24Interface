start_time = time()

using NaNStatistics
using LinearAlgebra, Random
using DelimitedFiles
using BSplineKit
using Printf

include("functions.jl");

box_length = 55.882675;
num_iter = 10 # 5000 ## pick a larger number

println("Reading pmf file and setting up optimization")
data = readdlm("freefile.txt", comments=true, comment_char='#');
xdata = data[:, 1];
ydata = data[:, 2];

fit_pmf = extrapolate(interpolate(xdata, ydata, BSplineOrder(6), Natural()), Flat());
fit_force = diff(fit_pmf, Derivative(1));

## parameters ( ref : https://pubs.acs.org/doi/pdf/10.1021/acs.jpca.2c03559 )
const kBT = 8.314 * 300.0; ## J/mol
const mass = 3.0 * 16e-3; ## kg/mol
const sig = sqrt(kBT / mass) * 1e-2; ## A/ps unit
const d_gas = 0.12 * 1e4;
const d_liquid = 1.76e-5 * 1e4; ## A^2/ps
const gamma_gas = kBT / mass / (d_gas) * 1e-4;
const gamma_liquid = kBT / mass / (d_liquid) * 1e-4; ## ps-1

function gamma_function(x, high_x=gamma_liquid, width_x=0.1, shift_x=20, power=1.5)
    result = @. (ifelse(x > 0, (tanh(width_x * (-x + shift_x)) / 2 + 0.5)^power * (high_x - gamma_gas) + gamma_gas, (tanh(width_x * (x + shift_x)) / 2 + 0.5)^power * (high_x - gamma_gas) + gamma_gas))
    return result
end;

# Setting up arrays
xmax = box_length;
const xlen = 401;
const vlen = 151;
xArr = range(-xmax, xmax, xlen);
vArr = range(-10.0 * sig, 10.0 * sig, vlen);

hx = abs(xArr[2] - xArr[1]);
hv = abs(vArr[2] - vArr[1]);
dt = 0.45 * min(hx, hv)^2;

x0 = gauss(xArr, 35.0, 2.0);
v0 = one_sided_gaussian(vArr, 0.0, sig, "L");

accuracy_order = 4;
if_periodic = false;
dx1 = 1.0 / hx * deriv1(xlen, accuracy_order, if_periodic);
dv1 = 1.0 / hv * deriv1(vlen, accuracy_order, false);
dv2 = 1.0 / hv / hv * deriv2(vlen, accuracy_order, false);

p0 = sparse(mychop(kron(x0, v0)));
pot_deriv = sparse(movmean(fit_force.(xArr), 10));
factor_force_over_mass = 4.184 * 1e-1
factor_kBT_over_mass = 1e-4

## time parameters
time_len = 30.0;
global p1 = p0;
dtt = dt / max(vlen, xlen);
println("Time step is ", dtt);
tlen = round(Integer, time_len / dtt);

global time_track = 0.0;
global p2 = zeros(length(p0));
global save_times = round.(Int, range(1.0, time_len, step=1.0) ./ dtt);

## get MD density data
println("Reading MD data files (", (time - start_time), " s)")
md_density = zeros(0)
for i in 1:round(Int,time_len)
    density_time = i
    pdata = readdlm("probability_data/final." * string(i) * ".txt", skipstart=1)
    # hdata = ifelse.(pdata[:, 2] .< box_length, pdata[:, 2], pdata[:, 2] .- 2.0 * box_length)
    hdata = pdata[:, 2]
    test_a = histcounts(hdata, xArr) / length(hdata)
    append!(md_density, test_a)
end
md_density = reshape(md_density, (xlen - 1, round(Int,time_len)+1));

md_velocity = zeros(0)
for i in 0:round(Int,time_len)
    density_time = i
    pdata = readdlm("probability_data/final." * string(i) * ".txt", skipstart=1)
    hdata = pdata[:, 3] .* 1e3 ## unit change to A/ps
    test_a = histcounts(hdata, vArr) / length(hdata)
    append!(md_velocity, test_a)
end
md_velocity = reshape(md_velocity, (vlen - 1, round(Int,time_len)+1));

if isfile("optimization_new.txt") == true
else
    file_v = open("optimization_new.txt","w")
    println(file_v, "#      Total Error (position)    "*"     param_1     "*"     param_2     "*"     param_3     ")
    close(file_v)
end

println("Normalizing density profiles from MD propagation")
for i in 1:round(Integer, total_time + 1)
    x = moving_average(xArr, 2)
    y = md_density[:, i]
    area = integrateArray(x, y)
    md_density[:,i] = md_density[:,i]/area
end
println("Normalizing density profiles from MD propagation")
for i in 1:round(Integer, total_time + 1)
    x = moving_average(vArr, 2)
    y = md_velocity[:, i]
    area = integrateArray(x, y)
    md_velocity[:,i] = md_velocity[:,i]/area
end

function optimize_loss(param_1::Float64, param_2::Float64, param_3::Float64) 
    global p1 = p0;    
    global time_track = 0.0;
    global p2 = zeros(length(p0));

    println("Input parameters are : ", param_1, " ", param_2, " ", param_3)
    gamma_array = sparse(gamma_function(xArr, gamma_liquid, param_1, param_2, param_3));

    function dpdt()
        term1 = -kron(dx1, spdiagm(vArr))
        term2 = kron(spdiagm(gamma_array), sparse(I, vlen, vlen)) + kron(spdiagm(gamma_array), spdiagm(vArr) * dv1)
        term3 = kron(spdiagm(pot_deriv) / mass * factor_force_over_mass, dv1)
        term4 = (kBT / mass) * factor_kBT_over_mass * kron(spdiagm(gamma_array), dv2)
        return @. (term1 + term2 + term3 + term4)
    end;
    dpdt_1 = dpdt();

    ## propagate
    global sum_error = 0
    # global sum_error_velocity = 0
    global kk = 0
    @fastmath @inbounds @simd for i = 1:tlen
        global time_track += dtt
        k1 = dpdt_1 * p1
        k2 = dpdt_1 * (p1 + k1 * 0.5 * dtt)
        k3 = dpdt_1 * (p1 + k2 * 0.5 * dtt)
        k4 = dpdt_1 * (p1 + k3 * dtt)
        @. global p2 = p1 + 1.0 / 6 * (k1 + k2 + k2 + k3 + k3 + k4) * dtt
        if i in save_times
            p0_dense = reshape(SparseArrays._SpecialArrays(p2), (vlen, xlen))
            x0_int = vec(sum(p0_dense, dims=1)) * hv
            v0_int = vec(sum(p0_dense, dims=2)) * hx

            area1 = hx * sum(x0_int)
            area2 = hv * sum(v0_int)

            x0_int = x0_int * area1
            v0_int = v0_int / area2
            
 	        kk += 1

            FP_x = moving_average(movmean(x0_int, 12), 2)
            FP_v = moving_average(movmean(v0_int, 3), 2)
            FP_x = FP_x / sum(FP_x)
            FP_v = FP_v / sum(FP_v)

            md_x = md_density[:,kk]
            md_x = md_x / sum(md_x)
            md_v = md_velocity[:,kk]
            md_v = md_v / sum(md_v)

	        local_error_position = sum(abs.( FP_x .- md_x ).^2)
            local_error_velocity = sum(abs.( FP_v .- md_v ).^2)
            global sum_error += local_error_position + local_error_velocity

            if isnan(area1) || isnan(area2) || isnan(sum_error)
                println("Check propagation parameters. Density is divergiving!")
                println("Local propagation error : ", local_error_position)
                println("sum(FP_x) is : ", sum(FP_x) )
                println("sum(md_x) is : ", sum(md_x) )
                println("sum(x0) : ", sum(x0_int))
                println("sum(v0) : ", sum(v0_int))
                println("Area 1 ", area1)
                println("Area2 ", area2)
                println("Exiting propagation loop...")
                exit()
            end;
            
	    println("Absolute error at t = ", kk - 1, " ps : ", local_error_position, " (in position) after ", (time()-start_time) , " seconds")
        end
        global p1 = p2
    end
    println("Total absolute error for the propagation : ", sum_error , " (in position) ")
    file_x = open("optimization_new.txt","a")
    @printf(file_x, "%.8f %.8f %.8f %.3f\n", sum_error, param_1, param_2, param_3)
    close(file_x)
    flush(stdout)
    return sum_error
end;


lower_bound = [0.15, 16.0, 0.8]
upper_bound = [0.8, 24.0, 3.0]

initial_guess = zeros(3)
for i=1:3
    initial_guess[i] = rand(lower_bound[i]:(upper_bound[i]-lower_bound[i])/100:upper_bound[i])
end

result = optimize(optimize_loss, lower_bound, upper_bound, initial_guess, SAMIN(nt=2, ns=2, rt=rand(range(0.85, 0.95, length=555)), verbosity=1), Optim.Options(show_trace=true, iterations=num_iter))
flush(stdout)

println("Done after ", (time() - start_time), " seconds")
