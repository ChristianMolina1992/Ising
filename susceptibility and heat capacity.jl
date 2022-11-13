using Random: default_rng, seed!
using RandomNumbers
using StaticArrays
using Plots
using BenchmarkTools

const β_crit = log(1+sqrt(2))/2
rand_ising2d(m, n=m) = rand(Int8[-1, 1], m, n)

seed!(4649)
s₀ = rand_ising2d(100);


function energy_density_ising2d(s)
    m, n = size(s)
    E = 0.0
    @inbounds begin
        for j in 1:n, i in 1:m-1
            E -= s[i,j]*s[i+1,j]
        end
        for j in 1:n
            E -= s[m,j]*s[1,j]
        end
        for j in 1:n-1, i in 1:m
            E -= s[i,j]*s[i,j+1]
        end
        for i in 1:m
            E -= s[i,m]*s[i,1]
        end
    end
    return E/(m*n)
end

function ising2d_ifelse!(s, β, niters, rng=default_rng())
    m, n = size(s)
    prob = @SVector [exp(-2*β*k) for k in -4:4]
    @fastmath @inbounds @simd for iter in 1:niters
        for j in 1:n 
            for i in 1:m
                let NN = s[ifelse(i == 1, m, i-1), j],  #ifelse(condition::Bool, x, y)Return x if condition is true, otherwise return y
                    SS = s[ifelse(i == m, 1, i+1), j],
                    WW = s[i, ifelse(j == 1, n, j-1)],
                    EE = s[i, ifelse(j == n, 1, j+1)],
                    CT = s[i, j]
                    k = CT * (NN + SS + WW + EE)
                    s[i,j] = ifelse(rand(rng) < prob[k+5], -CT, CT)
                end
            end
        end
    end
end

s = copy(s₀)
rng = default_rng()
ising2d_ifelse!(s, β_crit, 10, rng)

s = copy(s₀)
seed!(4649)
rng = default_rng()
@time ising2d_ifelse!(s, β_crit, 10^5, rng)
heatmap(s; size=(200, 200), axis=false, colorbar=false)
function magnetization_ising2d(s)
    m, n = size(s)
    return sum(s)/(m*n)
end
s = rand_ising2d(300)
magnetization_ising2d(s)    # Acá me tira error in method definition: function Ising2D.energy_density_ising2d must be explicitly imported to be extended



using Printf
using Plots
using Ising2D
using Statistics 


n_sweep = 100             # number of sweeps between sampling
n_therm = 1500           # number of sweeps to thermalize
n_data  = 100            # number of data samples per temperature
temps   = 4.0:-0.3:0.1   # temperatures to sample

function grafico()
    mt = Float64[] #Acá van a ir todos los puntos de la magnetización para cada T, los adjunto.
    et = Float64[]
    xt = Float64[]
    ct = Float64[]
    s = rand_ising2d(100) 
    for T in temps      #loop en la temperatura, termaliza
        ising2d_ifelse!(s,1/T,n_therm)
        m1=Float64[]
        e1=Float64[]
        
        for i=1:n_sweep
            ising2d_ifelse!(s,1/T,n_sweep)
            push!(m1,magnetization_ising2d(s))
            push!(e1,energy_density_ising2d(s))
                              
        end
        ma_ave = Statistics.mean(m1)  # Statistics.mean(m1)
        en_ave = Statistics.mean(e1)  # Statistics.mean(e1)
        susceptibility = Statistics.var(m1)/T
        heat_capacity = Statistics.var(e1)/T^2
                
        push!(mt,ma_ave)
        push!(et,en_ave)
        push!(xt,susceptibility)
        push!(ct,heat_capacity)
        
    end
     
    return collect(temps),xt
    
end
    
x= [2.4, 2.39, 2.38, 2.37, 2.36, 2.35, 2.34, 2.33, 2.32, 2.31, 2.3, 2.29, 2.28, 2.27, 2.26, 2.25, 2.24, 2.23, 2.22]          # Susceptibilidad vs Temperatura ( L = 100)
y=[0.006687215730639731, 0.0062480881433582686, 0.006858039215686276, 0.011237626184204918, 0.014085861605889406, 0.014233468471953576, 0.019851225546058886, 0.02618386209910261, 0.027769463504005578, 0.03953705417639599, 0.0547019501027668, 0.08113591836443033, 0.10204796306220097, 0.008315178162239131, 0.002949988198802182, 0.0018075724911335573, 0.0012135831529581527, 0.0015726794963083753, 0.0008292500846300846]
plot(x, y, xlabel = "T", ylabel = "Susceptibility", seriestype = :scatter)


x= [2.4, 2.39, 2.38, 2.37, 2.36, 2.35, 2.34, 2.33, 2.32, 2.31, 2.3, 2.29, 2.28, 2.27, 2.26, 2.25, 2.24, 2.23, 2.22]          # Susceptibilidad vs Temperatura ( L = 150)
y=[1.6497594144544542e-5, 1.9458437813335786e-5, 2.321310729030661e-5, 2.630397827041881e-5, 3.0314772495179537e-5, 3.637469098447108e-5, 4.788086355426073e-5, 5.493204652924119e-5, 8.26232911051771e-5, 0.00010825377322458884, 0.00016898758011766793, 0.000279193795335059, 0.0004713257057025261, 0.0006288446669553184, 0.0006529093955019728, 0.0009731466807087802, 5.280569342739598e-6, 3.4234149212863627e-6, 2.0104846222046647e-6]

x= [2.4, 2.39, 2.38, 2.37, 2.36, 2.35, 2.34, 2.33, 2.32, 2.31, 2.3, 2.29, 2.28, 2.27, 2.26, 2.25, 2.24, 2.23, 2.22]          # Susceptibilidad vs Temperatura ( L = 200)
y=[0.0017658561206860266, 0.0015255607889565106, 0.0025409535268652913, 0.0018249745262754122, 0.0028147799674713234, 0.005888550295615731, 0.006947181068268152, 0.006103225308340921, 0.008937428593151341, 0.017799347174231494, 0.010197619381971893, 0.029551756955912844, 0.04655927205653021, 0.0018472734449116715, 0.000601304505229284, 0.0006764143658810325, 0.000442268123647186, 0.00027481606864610236, 0.00016120554463554458]
plot(x, y, xlabel = "T", ylabel = "Susceptibility", seriestype = :scatter)


x= [2.4, 2.39, 2.38, 2.37, 2.36, 2.35, 2.34, 2.33, 2.32, 2.31, 2.3, 2.29, 2.28, 2.27, 2.26, 2.25, 2.24, 2.23, 2.22]          # Susceptibilidad vs Temperatura ( L = 250)
y=[0.0008078779741521885, 0.0013369618339038923, 0.0012737554233997115, 0.0016502981326541363, 0.001418019102986475, 0.0019198948354860947, 0.0036351933640680306, 0.00319202675569012, 0.004928095794647162, 0.00572773222151489, 0.011126992863440668, 0.010300404813914516, 0.00530367469540599, 0.003898781758151738, 0.0019012820486636278, 0.0011063715148598877, 0.00034663177027417026, 0.00015156555180069755, 0.00011816795860879055]
plot(x, y, xlabel = "T", ylabel = "Susceptibility", seriestype = :scatter)



x= [4.0, 3.8, 3.6, 3.4, 3.2, 3.0, 2.8, 2.6, 2.4, 2.2, 2.0, 1.8, 1.6, 1.4, 1.2, 1.0]          # Magnetización vs Temperatura ( L = 150)
y = [-0.0009795555555555556, 0.0018204444444444446, 0.00022488888888888897, -0.0029493333333333338, -0.002072, 7.199999999999985e-5, -0.0021182222222222223, 0.0008799999999999999, 0.020046222222222215, 0.7846426666666666, 0.912112888888889, 0.9565697777777779, 0.9794542222222222, 0.9912142222222222, 0.9970648888888887, 0.9993057777777778] 
plot(x, y, xlabel = "T", ylabel = "Magnetization", seriestype = :scatter)

x= [2.38, 2.36, 2.34, 2.32, 2.3, 2.28, 2.26, 2.24, 2.22, 2.2, 2.18, 2.16, 2.14, 2.12, 2.1, 2.08, 2.06, 2.04, 2.02]          # Magnetización vs Temperatura ( L = 150) problemas
y = [0.017060444444444448, -0.01163911111111111, 0.04601511111111111, 0.030666666666666672, -0.153192, 0.4903528888888888, 0.5695555555555556, 0.7115662222222221, 0.7594293333333334, 0.787656, 0.8084515555555555, 0.8251608888888887, 0.8439120000000001, 0.857104888888889, 0.8681217777777778, 0.8791022222222221, 0.8877120000000001, 0.8987697777777778, 0.9037022222222221] 
plot(x, y, xlabel = "T", ylabel = "Magnetization", seriestype = :scatter)


x= [4.0, 3.8, 3.6, 3.4, 3.2, 3.0, 2.8, 2.6, 2.4, 2.2, 2.0, 1.8, 1.6, 1.4, 1.2, 1.0]          # Energia vs Temperatura ( L = 150)
y = [-0.5567360000000001, -0.5952515555555555, -0.6375413333333334, -0.6856320000000001, -0.7466720000000001, -0.8166862222222222, -0.908312888888889, -1.028912, -1.202912, -1.547370666666667, -1.745416888888889, -1.860071111111111, -1.928981333333333, -1.9675715555555557, -1.9884284444444447, -1.9971448888888887] 
plot(x, y, xlabel = "T", ylabel = "Energy", seriestype = :scatter)



