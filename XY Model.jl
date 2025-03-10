using Random
using Statistics
using Plots

# Parámetros iniciales
T_values = Float64[]  # Temperaturas
mag_final = Float64[]  # Para almacenar la magnetización promedio en cada temperatura
susceptibilidad = Float64[]  # Para almacenar la susceptibilidad

pasos = 100000
m, n = 20, 20  # Tamaño de la matriz
s = rand(Uniform(0, 2pi), m, n)  # Inicialización de la matriz con ángulos aleatorios
N = m * n  # Número total de sitios

# Rango de temperaturas
temperaturas = 0.5:0.1:3.0  # Temperaturas de 1.0 a 4.0 con paso de 0.1

# Simulación de Monte Carlo para diferentes temperaturas
for T in temperaturas
    mag_temp = Float64[]  # Para almacenar la magnetización en cada iteración
    mag2_temp = Float64[]  # Para almacenar M^2 en cada iteración
    
    for iter in 1:pasos
        for j in 1:n 
            for i in 1:m
                # Calcular la energía actual del espín en (i, j)
                NN = s[ifelse(i == 1, m, i-1), j]
                SS = s[ifelse(i == m, 1, i+1), j]
                WW = s[i, ifelse(j == 1, n, j-1)]
                EE = s[i, ifelse(j == n, 1, j+1)]
                CT = s[i, j]
                E_before = -cos(CT - NN) - cos(CT - SS) - cos(CT - WW) - cos(CT - EE)

                # Proponer un nuevo ángulo para el espín en (i, j)
                new_angle = rand(Uniform(0, 2pi))
                E_after = -cos(new_angle - NN) - cos(new_angle - SS) - cos(new_angle - WW) - cos(new_angle - EE)

                # Calcular cambio de energía
                ΔE = E_after - E_before

                # Aceptar o rechazar el nuevo estado
                if ΔE < 0 || rand() < exp(-ΔE / T)
                    s[i, j] = new_angle
                end
            end
        end

        # Calcular la magnetización y M^2
        magnetizacion = sqrt(sum(cos.(s))^2 + sum(sin.(s))^2) / N
        push!(mag_temp, magnetizacion)
        push!(mag2_temp, magnetizacion^2)
    end

    # Promediar magnetización y M^2
    mean_m = mean(abs.(mag_temp))
    mean_m2 = mean(mag2_temp)
    
    # Calcular susceptibilidad
    chi = N * (mean_m2 - mean_m^2) / T

    # Almacenar los resultados
    push!(T_values, T)
    push!(mag_final, mean_m)
    push!(susceptibilidad, chi)
end

# Graficar magnetización vs temperatura
plot(T_values, mag_final, xlabel="Temperatura", ylabel="Magnetización", label="Magnetización", lw=2)

# Opcional: Graficar susceptibilidad vs temperatura
plot!(T_values, susceptibilidad, xlabel="Temperatura", ylabel="Susceptibilidad", label="Susceptibilidad", lw=2)










###################################

using Random
using Statistics
using Plots
using Distributions
using Base.Threads  # Para paralelización

# Parámetros iniciales
m, n = 20, 20  # Tamaño de la matriz
N = m * n  # Número total de sitios
pasos = 100000  # Número de pasos de Monte Carlo

# Rango de temperaturas
temperaturas = 0.5:0.1:3.0  # Temperaturas de 0.5 a 3.0 con paso de 0.1

# Arrays para almacenar resultados
T_values = Vector{Float64}(undef, length(temperaturas))  # Temperaturas
mag_final = Vector{Float64}(undef, length(temperaturas))  # Magnetización promedio
susceptibilidad = Vector{Float64}(undef, length(temperaturas))  # Susceptibilidad

# Simulación de Monte Carlo para diferentes temperaturas
@threads for (idx, T) in collect(enumerate(temperaturas))
    s = rand(Uniform(0, 2π), m, n)  # Inicialización de la matriz con ángulos aleatorios
    mag_temp = Vector{Float64}(undef, pasos)  # Para almacenar la magnetización en cada iteración
    mag2_temp = Vector{Float64}(undef, pasos)  # Para almacenar M^2 en cada iteración
    inv_T = 1 / T  # Precomputar 1 / T para evitar divisiones repetidas

    @inbounds for iter in 1:pasos
        for j in 1:n 
            for i in 1:m
                # Vecinos del espín en (i, j)
                NN = s[ifelse(i == 1, m, i-1), j]
                SS = s[ifelse(i == m, 1, i+1), j]
                WW = s[i, ifelse(j == 1, n, j-1)]
                EE = s[i, ifelse(j == n, 1, j+1)]
                CT = s[i, j]

                # Energía actual del espín en (i, j)
                E_before = -cos(CT - NN) - cos(CT - SS) - cos(CT - WW) - cos(CT - EE)

                # Proponer un nuevo ángulo para el espín en (i, j)
                new_angle = rand(Uniform(0, 2π))
                E_after = -cos(new_angle - NN) - cos(new_angle - SS) - cos(new_angle - WW) - cos(new_angle - EE)

                # Calcular cambio de energía
                ΔE = E_after - E_before

                # Aceptar o rechazar el nuevo estado
                if ΔE < 0 || rand() < exp(-ΔE * inv_T)
                    s[i, j] = new_angle
                end
            end
        end

        # Calcular la magnetización y M^2
        sum_cos = sum(cos, s)
        sum_sin = sum(sin, s)
        mag_temp[iter] = sqrt(sum_cos^2 + sum_sin^2) / N
        mag2_temp[iter] = mag_temp[iter]^2
    end

    # Promediar magnetización y M^2
    mean_m = mean(abs, mag_temp)
    mean_m2 = mean(mag2_temp)
    
    # Calcular susceptibilidad
    chi = N * (mean_m2 - mean_m^2) * inv_T

    # Almacenar los resultados
    T_values[idx] = T
    mag_final[idx] = mean_m
    susceptibilidad[idx] = chi
end

# Graficar magnetización vs temperatura
plot(T_values, mag_final, xlabel="Temperatura", ylabel="Magnetización", label="Magnetización", lw=2)

# Opcional: Graficar susceptibilidad vs temperatura
plot!(T_values, susceptibilidad, xlabel="Temperatura", ylabel="Susceptibilidad", label="Susceptibilidad", lw=2)
