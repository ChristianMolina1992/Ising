using Random, Plots, Statistics
gr()  # Activa el backend GR para Jupyter Notebook

# Parámetros del modelo
N = 32768          # Número de partículas
L = 256            # Tamaño del dominio
η = 0.5           # Nivel de ruido (0 = sin ruido, 1 = ruido máximo)
v = 0.5           # Velocidad de las partículas (debe ser escalar)
R = 1.0           # Radio de interacción
dt = 1.0          # Paso de tiempo
steps = 100       # Número de pasos de simulación

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

# Crear un objeto de animación
anim = Animation()

# Simulación con generación de imágenes
for step in 1:steps
    nuevas_θ = similar(θ)
    for i in 1:N
        vecinos = encontrar_vecinos(x, y, i, R)
        θ_prom = ángulo_promedio(θ, vecinos)
        nuevas_θ[i] = θ_prom + η * (rand() - 0.5) * 2π
    end

    θ .= nuevas_θ
    actualizar_posiciones!(x, y, θ, v, dt)  # Asegúrate que `v` es un escalar aquí

    # Guardar imagen cada 10 pasos
    if step % 10 == 0
        #println("Paso $step")

        # Preparar las coordenadas de las flechas
        u = cos.(θ)  # Componente x de la dirección de las flechas
        v_dir = sin.(θ)  # Componente y de la dirección de las flechas

        # Generar gráfico con flechas
        quiver(x, y, quiver=(u, v_dir), legend=false, xlim=(0, L), ylim=(0, L), 
               title="Paso $step", xlabel="Posición X", ylabel="Posición Y", 
               aspect_ratio=:equal, color=:blue)
        
        # Agregar el frame a la animación
        frame(anim)
    end
end

# Crear el GIF
gif(anim, "simulacion_flechas.gif", fps=1)  # Especifica la ruta y frames por segundo
  
