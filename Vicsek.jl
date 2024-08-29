##lo de los bloques

function asignar_a_celdas(x, y, L, r_cell)
    # Número de celdas por lado
    n_celdas = ceil(Int, L / r_cell)
    
    # Inicializar las celdas vacías
    celdas = [Int[] for _ in 1:n_celdas, _ in 1:n_celdas]
    
    # Asignar cada partícula a una celda
    for i in 1:length(x)
        ix = cld(x[i], r_cell)
        iy = cld(y[i], r_cell)
        push!(celdas[ix, iy], i)
    end
    
    return celdas, n_celdas
end

function encontrar_vecinos_por_celdas(x, y, i, R, r_cell, celdas, n_celdas)
    vecinos = []
    
    # Determinar la celda de la partícula i
    ix = cld(x[i], r_cell)
    iy = cld(y[i], r_cell)
    
    # Recorrer las celdas vecinas
    for dx in -1:1
        for dy in -1:1
            # Calcular las celdas vecinas considerando periodicidad
            ix_vecina = mod1(ix + dx, n_celdas)
            iy_vecina = mod1(iy + dy, n_celdas)
            
            # Revisar todas las partículas en la celda vecina
            for j in celdas[ix_vecina, iy_vecina]
                if i != j
                    dx = abs(x[i] - x[j])
                    dy = abs(y[i] - y[j])
                    dx = min(dx, L - dx)
                    dy = min(dy, L - dy)
                    if dx^2 + dy^2 < R^2
                        push!(vecinos, j)
                    end
                end
            end
        end
    end
    
    return vecinos
end


##Lo otro 

using Random, Plots, Statistics
gr()  # Activa el backend GR para Jupyter Notebook

# Parámetros del modelo
N = 1000          # Número de partículas
L = 10.0          # Tamaño del dominio
η = 0.1           # Nivel de ruido (0 = sin ruido, 1 = ruido máximo)
v = 0.03          # Velocidad de las partículas
R = 1.0           # Radio de interacción
dt = 1.0          # Paso de tiempo
steps = 1000      # Número de pasos de simulación

# Inicialización de las partículas
x = L .* rand(N)  # Posiciones x
y = L .* rand(N)  # Posiciones y
θ = 2π .* rand(N) # Ángulos de dirección

# Función para actualizar las posiciones
function actualizar_posiciones!(x, y, θ, v, dt)
    x .+= v * cos.(θ) * dt
    y .+= v * sin.(θ) * dt

    # Condiciones de frontera periódicas
    x .= mod.(x, L)
    y .= mod.(y, L)
end

# Función para calcular el ángulo promedio de los vecinos
function ángulo_promedio(θ, vecinos)
    sin_prom = mean(sin.(θ[vecinos]))
    cos_prom = mean(cos.(θ[vecinos]))
    return atan(sin_prom, cos_prom)
end

# Función para encontrar vecinos
function encontrar_vecinos(x, y, i, R)
    vecinos = []
    for j in 1:N
        if i != j
            dx = abs(x[i] - x[j])
            dy = abs(y[i] - y[j])
            # Considera condiciones de frontera periódicas
            dx = min(dx, L - dx)
            dy = min(dy, L - dy)
            if dx^2 + dy^2 < R^2
                push!(vecinos, j)
            end
        end
    end
    return vecinos
end

# Simulación
for step in 1:steps
    nuevas_θ = similar(θ)
    for i in 1:N
        vecinos = encontrar_vecinos(x, y, i, R)
        θ_prom = ángulo_promedio(θ, vecinos)
        nuevas_θ[i] = θ_prom + η * (rand() - 0.5) * 2π
    end

    θ .= nuevas_θ
    actualizar_posiciones!(x, y, θ, v, dt)

    if step % 100 == 0
        p = plot(x, y, seriestype = :scatter, markersize = 2, aspect_ratio = 1, legend = false)
        title!(p, "Paso $step")
        xlims!(0, L)
        ylims!(0, L)
        display(p)
    end
end


