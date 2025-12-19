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




###WOLFF###

using Random
using Statistics
using Plots

# ============================================================
# 1) Inicialización del lattice
# ============================================================
function inicializar_lattice(L)
    # ángulos uniformes en [0, 2π)
    return 2π .* rand(L, L)
end

# ============================================================
# 2) Magnetización vectorial
# ============================================================
function magnetizacion(theta)
    mx = mean(cos.(theta))
    my = mean(sin.(theta))
    return sqrt(mx^2 + my^2)
end

# ============================================================
# 3) Wolff step CORRECTO para modelo XY
# ============================================================
function wolff_step!(theta, β)
    L = size(theta, 1)

    # --------------------------------------------------------
    # 1) Dirección aleatoria para el hiperplano
    # --------------------------------------------------------
    φ = 2π * rand()
    nx, ny = cos(φ), sin(φ)

    # Proyección del spin sobre la dirección n
    proj(i, j) = nx * cos(theta[i, j]) + ny * sin(theta[i, j])

    # --------------------------------------------------------
    # 2) Elegir semilla del cluster
    # --------------------------------------------------------
    i0 = rand(1:L)
    j0 = rand(1:L)

    in_cluster = falses(L, L)
    cluster    = [(i0, j0)]
    stack      = [(i0, j0)]

    in_cluster[i0, j0] = true
    S0 = proj(i0, j0)

    # --------------------------------------------------------
    # 3) Construcción del cluster
    # --------------------------------------------------------
    while !isempty(stack)
        i, j = pop!(stack)
        Si = proj(i, j)

        vecinos = (
            (mod1(i + 1, L), j),
            (mod1(i - 1, L), j),
            (i, mod1(j + 1, L)),
            (i, mod1(j - 1, L))
        )

        for (a, b) in vecinos
            if !in_cluster[a, b]
                Sn = proj(a, b)

                # Solo enlaces ferromagnéticos proyectados
                if Si * Sn > 0
                    p_add = 1 - exp(-2β * Si * Sn)

                    if rand() < p_add
                        in_cluster[a, b] = true
                        push!(cluster, (a, b))
                        push!(stack, (a, b))
                    end
                end
            end
        end
    end

    # --------------------------------------------------------
    # 4) Reflección correcta de todo el cluster
    # --------------------------------------------------------
    for (i, j) in cluster
        sx = cos(theta[i, j])
        sy = sin(theta[i, j])

        dot = sx * nx + sy * ny   # S · n

        sx2 = sx - 2 * dot * nx
        sy2 = sy - 2 * dot * ny

        theta[i, j] = atan(sy2, sx2)
    end

    return nothing
end

function simular_wolff(;L, T,pasos_equil, pasos_med)

    mags = Float64[]
    chis = Float64[]

    β = 1 / T
    theta = inicializar_lattice(L)

    # ---------------- Equilibración ----------------
    for _ in 1:pasos_equil
        wolff_step!(theta, β)
    end

    # ---------------- Medición ---------------------
    medidas  = Float64[]
    medidas2 = Float64[]

    for _ in 1:pasos_med
        wolff_step!(theta, β)
        m = magnetizacion(theta)
        push!(medidas, m)
        push!(medidas2, m^2)
    end

    m1 = mean(medidas)
    m2 = mean(medidas2)

    chi = L^2 * (m2 - m1^2) / T

    push!(mags, m1)
    push!(chis, chi)

    println("Wolff: L = $L  T = $T → m = $m1   χ = $chi")

    return m1, chi, theta
end

using DelimitedFiles


L = 100
T = 0.89290

m, chi, theta = simular_wolff(L = L,T = T,pasos_equil = 100000,pasos_med  = 20000000)

mkpath("/mnt/d/Doctorado/xymodel/L=$(L)")
writedlm("/mnt/d/Doctorado/xymodel/L=$(L)/theta_final.txt", theta)

using DelimitedFiles

# ============================================================
# 5) Generar N snapshots, separados por Δ pasos Wolff
# ============================================================
function generar_snapshots_wolff(;
    L,
    T,
    n_archivos,          # cantidad de archivos que querés
    pasos_entre,         # pasos Wolff entre archivos (ej: 10000)
    carpeta_base = "/mnt/d/Doctorado/xymodel"
)

    β = 1 / T

    # --------------------------------------------------------
    # 1) Cargar configuración inicial
    # --------------------------------------------------------
    archivo_theta = "$(carpeta_base)/L=$(L)/theta_final.txt"
    theta = readdlm(archivo_theta)

    println("Configuración inicial cargada:")
    println(archivo_theta)

    pasos_acumulados = 0

    # --------------------------------------------------------
    # 2) Loop principal
    # --------------------------------------------------------
    for k in 1:n_archivos

        # correr Wolff pasos_entre veces
        for _ in 1:pasos_entre
            wolff_step!(theta, β)
        end

        pasos_acumulados += pasos_entre

        nombre = "$(carpeta_base)/L=$(L)/XYconfig_L$(L)_pasos$(pasos_acumulados).txt"
        writedlm(nombre, theta)

    end

    return theta
end
L = 100
T = 0.89290

theta_final = generar_snapshots_wolff(
    L = L,
    T = T,
    n_archivos  = 100000,
    pasos_entre = 10000
)
