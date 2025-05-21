using DelimitedFiles
using Plots

# Parámetros críticos del modelo de Ising 2D
Tc = Ising_SQ_critical_temperature  # Temperatura crítica
γ = 7/4                             # Exponente de susceptibilidad
ν = 1                               # Exponente de correlación

# Función para leer datos desde archivo
function leer_datos(L)
    ruta = "/media/cmolina/Christian1/Doctorado/L=$(L)/magnetizacion_susceptibilidad_L$(L).txt"
    datos = readdlm(ruta, skipstart=1)  # Saltar encabezado

    T = datos[:, 1]
    M = datos[:, 2]
    χ = datos[:, 3]

    return T, M, χ
end

# Lista de tamaños de red
Ls = [20, 50, 100]

# Inicializar gráfico
plot(
    title = "Colapso de susceptibilidad Ising 2D",
    xlabel = "(T - Tc)/Tc × L^{1/ν}",
    ylabel = "χ × L^{-γ/ν}",
    legend = :topleft,
    lw = 2,
    grid = true
)

# Bucle sobre cada tamaño L
for L in Ls
    T, _, χ = leer_datos(L)

    # Variable de escala reducida
    x = ((T .- Tc) ./ Tc) .* L^(1/ν)
    y = χ .* L^(-γ/ν)

    plot!(x, y, label = "L = $L", marker = :auto)
end

# Mostrar gráfico
display(current())

# Guardar (opcional)
savefig("colapso_susceptibilidad_reducida_ising2D.png")




using DelimitedFiles
using Plots

# Parámetros críticos del modelo de Ising 2D
Tc = Ising_SQ_critical_temperature  # Temperatura crítica
γ = 7/4                             # Exponente de susceptibilidad
ν = 1                               # Exponente de correlación

# Función para leer datos desde archivo
function leer_datos(L)
    ruta = "/media/cmolina/Christian1/Doctorado/L=$(L)/magnetizacion_susceptibilidad_L$(L).txt"
    datos = readdlm(ruta, skipstart=1)  # Saltar encabezado

    T = datos[:, 1]
    M = datos[:, 2]
    χ = datos[:, 3]

    return T, M, χ
end

# Lista de tamaños de red
Ls = [20, 50, 100]

# Inicializar gráfico
plot(
    title = "Colapso de susceptibilidad Ising 2D",
    xlabel = "(T - Tc)/Tc × L^{1/ν}",
    ylabel = "χ × L^{-γ/ν}",
    legend = :topleft,
    lw = 2,
    grid = true
)

# Bucle sobre cada tamaño L
for L in Ls
    T, _, χ = leer_datos(L)

    # Variable de escala reducida
    x = ((T .- Tc) ./ Tc) .* L^(1/ν)
    y = χ .* L^(-γ/ν)

    plot!(x, y, label = "L = $L", marker = :auto)
end

# Mostrar gráfico
display(current())

# Guardar (opcional)
savefig("colapso_susceptibilidad_reducida_ising2D.png")




###########ESto es la primera parte
using Plots
using DelimitedFiles
using Statistics

# Temperatura crítica (Ising 2D)
Tc = 2.0 / log(1.0 + sqrt(2.0))  

# Función para cargar datos
function cargar_datos(L)
    ruta = "/media/cmolina/Christian1/Doctorado/Prueba/L=$(L)/magnetizacion_susceptibilidad_L$(L).txt"
    datos = readdlm(ruta, skipstart=1)
    return datos[:, 1], datos[:, 2], datos[:, 3]  # T, mag, sus
end

# Cargar datos
Ls = [20, 50, 100, 200]
datos = Dict{Int, Any}()

for L in Ls
    T, mag, sus = cargar_datos(L)
    datos[L] = (T=T, mag=mag, sus=sus)
end

# Parámetros de escalamiento
β = 1/8    # Exponente de magnetización
γ = 7/4    # Exponente de susceptibilidad
ν = 1      # Exponente de longitud de correlación

# --- Gráfico de magnetización escalada ---
p1 = plot(xlabel="T/Tc", ylabel="M · L^{β/ν}", title="Magnetización escalada")
for L in Ls
    T = datos[L].T
    scaled_mag = datos[L].mag .* (L^(β/ν))
    plot!((T.-Tc)/Tc, scaled_mag, label="L=$L", lw=2)
end
vline!([1.0], label="Tc", color=:black, ls=:dash)

# --- Gráfico de susceptibilidad escalada ---
p2 = plot(
    xlabel="(T - Tc) · L^{1/ν} / Tc", 
    ylabel="χ / L^{γ/ν}", 
    title="Susceptibilidad escalada (colapso finito)"
)
for L in Ls
    T = datos[L].T
    scaled_T = (T .- Tc) .* (L^(1/ν)) ./ Tc  # Escalamiento correcto
    scaled_sus = datos[L].sus ./ (L^(γ/ν))
    plot!(scaled_T, scaled_sus, label="L=$L", lw=2, marker=:circle, ms=3)
end
vline!([0.0], label="Tc", color=:black, ls=:dash)

# Mostrar gráficos
plot(p1, p2, layout=(2,1), size=(900, 1200), legend=:topright)


using Plots
using DelimitedFiles
using Statistics
using LinearAlgebra  # Para ajuste lineal (opcional)

# Cargar datos (asegúrate de que esta parte ya está ejecutada)
# datos[L] = (T=..., mag=..., sus=...)

# --- Función para encontrar T_c(L) a partir del máximo de susceptibilidad ---
function encontrar_Tc_L(datos, Ls)
    Tc_L = Dict{Int, Float64}()
    for L in Ls
        T = datos[L].T
        sus = datos[L].sus
        idx_max = argmax(sus)  # Índice del máximo de susceptibilidad
        Tc_L[L] = T[idx_max]   # Temperatura crítica para este L
    end
    return Tc_L
end

# Calcular Tc(L) para cada L
Ls = [20, 50, 100, 200]
Tc_L = encontrar_Tc_L(datos, Ls)

# --- Gráfico de susceptibilidad con máximos marcados ---
p = plot(xlabel="T", ylabel="Susceptibilidad (χ)", title="Susceptibilidad por tamaño L", legend=:topright)

for L in Ls
    T = datos[L].T
    sus = datos[L].sus
    plot!(T, sus, label="L=$L", lw=2)
    
    # Marcar el punto de Tc(L)
    scatter!([Tc_L[L]], [maximum(sus)], color=:red, marker=:star, label="Tc(L=$L) = $(round(Tc_L[L], digits=4))")
end

# Añadir línea vertical en Tc exacta (opcional)
Tc = 2.0 / log(1.0 + sqrt(2.0))
vline!([Tc], label="Tc exacta = $(round(Tc, digits=4))", color=:black, linestyle=:dash)

# Mostrar gráfico
display(p)

# --- Resultados numéricos ---
println("\nTemperaturas críticas estimadas:")
for L in Ls
    println("L = $L: Tc(L) = $(round(Tc_L[L], digits=5))")
end

# --- Extrapolación a L → ∞ (opcional) ---
if length(Ls) >= 2
    L_inv = 1 ./ Ls
    Tc_values = [Tc_L[L] for L in Ls]
    
    # Ajuste lineal: Tc(L) ≈ Tc_inf + a*(1/L)
    A = hcat(L_inv, ones(length(L_inv)))
    coeffs = A \ Tc_values
    Tc_inf = coeffs[2]  # Tc en L → ∞
    
    println("\nExtrapolación a L → ∞: Tc = $(round(Tc_inf, digits=5))")
end


using JLD
using LatticeModels
using DelimitedFiles
using Statistics
using BioStatPhys
using Ising2D

# Temperatura crítica
Tc = Ising_SQ_critical_temperature

# Rango de temperaturas (proporcional a Tc)
Temperaturas_modulo = collect(0.85:0.01:1.15) * Tc

# Función para generar datos a una temperatura T y tamaño L
function config(; T, L, nterm, n)
    IS = Ising(SQLattice_periodic, L, L, ordered=true)

    set_temperature!(IS, T)
    Metropolis!(IS, steps=nterm)
    _, m = Metropolis!(IS, steps=n, save_interval=10)

    mag = Statistics.mean(abs.(m))
    sus = (L * L * Statistics.var(abs.(m))) / T

    # === Guardar configuración final ===
    carpeta = "/media/cmolina/Christian1/Doctorado/Prueba/L=$(L)/"
    mkpath(carpeta)

    nombrebase = "SQconfig_L$(L)_Temp$(round(T, digits=5))_ndata$(n)"
    archivo_txt = carpeta * nombrebase * ".txt"
    archivo_jld = carpeta * nombrebase * ".jld"

    writedlm(archivo_txt, IS.σ)   # Guardar en texto plano
    @save archivo_jld IS.σ        # Guardar en formato binario .jld

    return mag, sus
end

# Función para calcular magnetización y susceptibilidad para un tamaño L
function sus_mag(; L, nterm=700000, n=700000)
    magnetizacion = Float64[]
    susceptibilidad = Float64[]

    for T in Temperaturas_modulo
        mag, sus = config(T=T, L=L, nterm=nterm, n=n)
        push!(magnetizacion, mag)
        push!(susceptibilidad, sus)
    end
    return magnetizacion, susceptibilidad
end

# Función para guardar los resultados en archivo
function guardar_datos(L, temperaturas, magnetizacion, susceptibilidad)
    datos = hcat(temperaturas, magnetizacion, susceptibilidad)
    ruta = "/media/cmolina/Christian1/Doctorado/Prueba/L=$(L)/magnetizacion_susceptibilidad_L$(L).txt"
    mkpath(dirname(ruta))

    open(ruta, "w") do io
        write(io, "Temperatura Magnetizacion Susceptibilidad\n")
        writedlm(io, datos)
    end
end

# === Bucle principal para todos los tamaños L deseados ===
for L in [20,50,100,200]
    println("Procesando L = $L")
    magnetizacion, susceptibilidad = sus_mag(L=L)
    guardar_datos(L, Temperaturas_modulo, magnetizacion, susceptibilidad)
end


