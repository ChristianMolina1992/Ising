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
