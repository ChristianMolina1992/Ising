using Random, Plots, Statistics
gr()  # Activa el backend GR para Jupyter Notebook

# Parámetros del modelo
N = 50          # Número de partículas
L = 20            # Tamaño del dominio
η = 0           # Nivel de ruido (0 = sin ruido, 1 = ruido máximo)
v = 1.0           # Velocidad de las partículas (debe ser escalar)
R = 1.0           # Radio de interacción
dt = 1.0          # Paso de tiempo
steps = 2000       # Número de pasos de simulación

# Inicialización de las partículas
x = L .* rand(N)  # Posiciones x
y = L .* rand(N)  # Posiciones y
θ = 2π .* rand(N) # Ángulos de dirección

# Función para actualizar las posiciones
function actualizar_posiciones!(x, y, θ, v, dt)
    x .+= v * cos.(θ) * dt
    y .+= v * sin.(θ) * dt

    # Aplicar condiciones de frontera periódicas
    x .= (x .% L) .+ (x .< 0) * L
    y .= (y .% L) .+ (y .< 0) * L
end

# Función para calcular el ángulo promedio de los vecinos
function ángulo_promedio(θ, vecinos)
    if length(vecinos) == 0
        return θ[1]  # Retornar el ángulo actual si no hay vecinos
    end
    sin_prom = mean(sin.(θ[vecinos]))
    cos_prom = mean(cos.(θ[vecinos]))
    return atan(sin_prom, cos_prom)
end

function asignar_celdas(x, y, L, R)
    # Crear un diccionario para almacenar las partículas en cada celda
    celdas = Dict{Tuple{Int, Int}, Vector{Int}}()
    for i in 1:length(x)
        # Calcular la celda en la que está la partícula i
        celda_x = Int(floor(x[i] / R))
        celda_y = Int(floor(y[i] / R))
        celda = (celda_x, celda_y)
        # Usar `haskey` para verificar si la celda ya está en el diccionario
        if !haskey(celdas, celda)
            celdas[celda] = []
        end
        push!(celdas[celda], i)
    end
    return celdas
end

function encontrar_vecinos_celdas(x, y, i, R, L, celdas)
    vecinos = []
    # Calcular la celda en la que está la partícula i
    celda_x = Int(floor(x[i] / R))
    celda_y = Int(floor(y[i] / R))
    # Revisar partículas en la misma celda y en las celdas vecinas
    for dx in -1:1, dy in -1:1
        celda_vecina = (mod(celda_x + dx, L ÷ R), mod(celda_y + dy, L ÷ R))
        if haskey(celdas, celda_vecina)  # Cambiado a haskey
            for j in celdas[celda_vecina]
                if i != j
                    # Calcular la distancia entre partículas con condiciones periódicas
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


# Medir el tiempo de ejecución
@time begin
    # Crear un objeto de animación
    anim = Animation()

    for step in 1:steps
        # Asignar partículas a las celdas en cada paso
        celdas = asignar_celdas(x, y, L, R)

        nuevas_θ = similar(θ)
        for i in 1:N
            vecinos = encontrar_vecinos_celdas(x, y, i, R, L, celdas)
            θ_prom = ángulo_promedio(θ, vecinos)
            nuevas_θ[i] = θ_prom + η * (rand() - 0.5) * 2π
        end

        θ .= nuevas_θ
        actualizar_posiciones!(x, y, θ, v, dt)

        # Guardar imagen cada 10 pasos
        if step % 10 == 0
            println("Paso $step")

            # Generar gráfico solo con posiciones
            scatter(x, y, legend=false, xlim=(0, L), ylim=(0, L), 
                    title="Paso $step", xlabel="Posición X", ylabel="Posición Y", 
                    color=:blue, markersize=5)
            
            # Agregar el frame a la animación
            frame(anim)
        end
    end

    # Crear el GIF
    gif(anim, "simulacion_posiciones.gif", fps=0.8)  # Especifica la ruta y frames por segundo
end



using Random, Plots, Statistics
gr()  # Activa el backend GR para Jupyter Notebook

# Parámetros del modelo
N = 32768          # Número de partículas
L = 256          # Tamaño del dominio
v = 0.3         # Velocidad de las partículas (debe ser escalar)
R = 1.0         # Radio de interacción
dt = 1.0        # Paso de tiempo
steps = 2000    # Número de pasos de simulación
η_values = 0:0.1:1  # Valores de eta (ruido)

# Función para calcular el parámetro de orden
function calcular_parametro_orden(θ)
    # Suma de los vectores velocidad
    v_x = sum(cos.(θ))
    v_y = sum(sin.(θ))
    # Magnitud de la suma dividida por el número de partículas
    return sqrt(v_x^2 + v_y^2) / N
end

# Inicialización de las partículas
x = L .* rand(N)  # Posiciones x
y = L .* rand(N)  # Posiciones y

# Contenedor para los resultados del parámetro de orden
parametro_orden_vs_eta = []

for η in η_values
    θ = 2π .* rand(N) # Ángulos de dirección para cada nuevo η

    for step in 1:steps
        celdas = asignar_celdas(x, y, L, R)

        nuevas_θ = similar(θ)
        for i in 1:N
            vecinos = encontrar_vecinos_celdas(x, y, i, R, L, celdas)
            θ_prom = ángulo_promedio(θ, vecinos)
            nuevas_θ[i] = θ_prom + η * (rand() - 0.5) * 2π
        end

        θ .= nuevas_θ
        actualizar_posiciones!(x, y, θ, v, dt)
    end

    # Calcular el parámetro de orden para el valor de η actual
    parametro_orden = calcular_parametro_orden(θ)
    push!(parametro_orden_vs_eta, parametro_orden)
end

# Graficar el parámetro de orden en función de η
plot(η_values, parametro_orden_vs_eta, xlabel="η (Ruido)", ylabel="Parámetro de Orden (velocidades)", 
     label="L=64", lw=2, legend=:top)
  
