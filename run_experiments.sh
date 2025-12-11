#!/bin/bash
# ==========================
# run_experiments.sh
# Ejecuta instancias SECUENCIALMENTE para FILO2
# Organiza resultados por configuración
# caffeinate -i -s ./run_experiments.sh
# ==========================
 
# CONFIGURACIÓN - Modifica según tus necesidades
# ==================================================

# DATASETS=("X" "B" "I")
# SEEDS=(0 1 2 3 4 5 6 7 8 9)
# TIME_LIMIT=300

DATASETS=("X" "B")  # Datasets a procesar
SEEDS=(0 1 2 3 4)       # Seeds a ejecutar
TIME_LIMIT=180          # Tiempo límite en segundos (5 min)

# Nombre de esta configuración experimental (para organizar resultados)
EXPERIMENT_NAME="XB-180"

# Parámetros adicionales del solver (opcional)
SOLVER_PARAMS="--time-limit $TIME_LIMIT --save-trajectory"
# Ejemplos con más parámetros:
# SOLVER_PARAMS="--time-limit $TIME_LIMIT --save-trajectory --granular-neighbors 30"
# SOLVER_PARAMS="--unit-demands --time-limit $TIME_LIMIT --save-trajectory"

# ==================================================

# Colores para output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Directorio base de resultados
RESULTS_BASE="results/$EXPERIMENT_NAME"
RESULTS_DIR="$RESULTS_BASE/outputs"
LOGS_DIR="$RESULTS_BASE/logs"
TRAJECTORIES_DIR="$RESULTS_BASE/trajectories"

# Crear estructura de directorios
mkdir -p "$RESULTS_DIR"
mkdir -p "$LOGS_DIR"
mkdir -p "$TRAJECTORIES_DIR"

# Archivo de resumen
SUMMARY_FILE="$RESULTS_BASE/summary.csv"
echo "instance,seed,cost,time,routes,trajectory_file" > "$SUMMARY_FILE"

# Contador de progreso
total_instances=0
completed_instances=0

# Primero contar cuántas instancias hay
echo -e "${BLUE}Contando instancias...${NC}"
for dataset in "${DATASETS[@]}"; do
    if [ "$dataset" = "B" ]; then
        EXT="txt"
    else
        EXT="vrp"
    fi
    
    for inst in ../instances/$dataset/*.$EXT; do
        if [ -f "$inst" ]; then
            total_instances=$((total_instances + ${#SEEDS[@]}))
        fi
    done
done

echo -e "${GREEN}Total de ejecuciones a realizar: $total_instances${NC}"
echo -e "${YELLOW}Configuración: $EXPERIMENT_NAME${NC}"
echo -e "${YELLOW}Time limit: ${TIME_LIMIT}s por instancia${NC}"
echo ""

# Timestamp de inicio
start_time=$(date +%s)
echo "Inicio: $(date)" | tee "$RESULTS_BASE/experiment_info.txt"
echo "Configuración: $EXPERIMENT_NAME" >> "$RESULTS_BASE/experiment_info.txt"
echo "Time limit: ${TIME_LIMIT}s" >> "$RESULTS_BASE/experiment_info.txt"
echo "Seeds: ${SEEDS[*]}" >> "$RESULTS_BASE/experiment_info.txt"
echo "Parámetros: $SOLVER_PARAMS" >> "$RESULTS_BASE/experiment_info.txt"
echo "" >> "$RESULTS_BASE/experiment_info.txt"

# Función para estimar tiempo restante
estimate_remaining_time() {
    current_time=$(date +%s)
    elapsed=$((current_time - start_time))
    
    if [ $completed_instances -gt 0 ]; then
        avg_time=$((elapsed / completed_instances))
        remaining=$((total_instances - completed_instances))
        eta=$((remaining * avg_time))
        
        # Convertir a formato legible
        hours=$((eta / 3600))
        minutes=$(((eta % 3600) / 60))
        seconds=$((eta % 60))
        
        echo "${hours}h ${minutes}m ${seconds}s"
    else
        echo "Calculando..."
    fi
}

# Procesar cada dataset
for dataset in "${DATASETS[@]}"; do
    echo -e "\n${BLUE}========================================${NC}"
    echo -e "${BLUE}Procesando dataset: $dataset${NC}"
    echo -e "${BLUE}========================================${NC}\n"
    
    # Definir extensión
    if [ "$dataset" = "B" ]; then
        EXT="txt"
    else
        EXT="vrp"
    fi
    
    # Procesar cada instancia
    for inst in ../instances/$dataset/*.$EXT; do
        if [ ! -f "$inst" ]; then
            continue
        fi
        
        # Nombre base de la instancia
        inst_name=$(basename "$inst")
        inst_name_noext="${inst_name%.*}"
        
        echo -e "${GREEN}Instancia: $inst_name_noext${NC}"
        
        # Ejecutar para cada seed
        for seed in "${SEEDS[@]}"; do
            completed_instances=$((completed_instances + 1))
            
            # Nombres de archivos (ajustado al formato real del solver)
            log_file="$LOGS_DIR/${inst_name_noext}_seed${seed}.log"
            
            # El solver genera nombres con el formato: basename_seed-N.ext
            # donde basename incluye la extensión original (.vrp, .txt)
            out_file_src="${inst_name}_seed-${seed}.out"
            sol_file_src="${inst_name}_seed-${seed}.vrp.sol"
            traj_file_src="${inst_name}_seed-${seed}.trajectory"
            
            # Progress bar
            progress=$((completed_instances * 100 / total_instances))
            eta=$(estimate_remaining_time)
            
            echo -ne "  [${completed_instances}/${total_instances}] Seed $seed - "
            echo -ne "Progress: ${progress}% - ETA: ${eta}   \r"
            
            # Ejecutar el solver (con --outpath para controlar donde se generan los archivos)
            ./filo2 "$inst" --seed $seed --outpath "$RESULTS_DIR/" $SOLVER_PARAMS > "$log_file" 2>&1
            
            # Los archivos ahora se generan directamente en RESULTS_DIR
            out_file="$RESULTS_DIR/$out_file_src"
            sol_file="$RESULTS_DIR/$sol_file_src"
            traj_file="$RESULTS_DIR/$traj_file_src"
            
            # Verificar que se generaron los archivos
            if [ -f "$out_file" ]; then
                # Leer resultados del .out
                cost=$(awk '{print $1}' "$out_file")
                time=$(awk '{print $2}' "$out_file")
                
                # Leer número de rutas del log
                # Linux
                # routes=$(grep -oP 'n\. routes = \K\d+' "$log_file" | tail -1)
                # macOS
                routes=$(grep 'n\. routes =' "$log_file" | tail -1 | sed 's/.*n\. routes = \([0-9]*\).*/\1/')
                
                # Mover trajectory a su carpeta si existe
                if [ -f "$traj_file" ]; then
                    mv "$traj_file" "$TRAJECTORIES_DIR/"
                    traj_path="$TRAJECTORIES_DIR/$traj_file_src"
                else
                    traj_path="NA"
                fi
                
                # Agregar al resumen
                echo "$inst_name_noext,$seed,$cost,$time,$routes,$traj_path" >> "$SUMMARY_FILE"
                
                echo -ne "  [${completed_instances}/${total_instances}] Seed $seed - "
                echo -e "${GREEN}✓${NC} Cost: $cost Time: ${time}s Routes: $routes"
            else
                echo -ne "  [${completed_instances}/${total_instances}] Seed $seed - "
                echo -e "${RED}✗ ERROR${NC} (ver log: $log_file)"
                echo "$inst_name_noext,$seed,ERROR,ERROR,ERROR,NA" >> "$SUMMARY_FILE"
                
                # Debug: mostrar qué archivos se generaron
                echo "    Debug: Archivos generados:"
                ls -la "$RESULTS_DIR/" | grep "${inst_name_noext}" | tail -3
            fi
        done
        
        echo ""
    done
done

# Timestamp de fin
end_time=$(date +%s)
total_time=$((end_time - start_time))
hours=$((total_time / 3600))
minutes=$(((total_time % 3600) / 60))
seconds=$((total_time % 60))

echo -e "\n${GREEN}========================================${NC}"
echo -e "${GREEN}EXPERIMENTO COMPLETADO${NC}"
echo -e "${GREEN}========================================${NC}"
echo -e "Tiempo total: ${hours}h ${minutes}m ${seconds}s"
echo -e "Instancias procesadas: ${completed_instances}"
echo -e "\nResultados guardados en: ${BLUE}$RESULTS_BASE${NC}"
echo -e "  - Outputs: ${BLUE}$RESULTS_DIR${NC}"
echo -e "  - Logs: ${BLUE}$LOGS_DIR${NC}"
echo -e "  - Trajectories: ${BLUE}$TRAJECTORIES_DIR${NC}"
echo -e "  - Resumen: ${BLUE}$SUMMARY_FILE${NC}"

echo "" >> "$RESULTS_BASE/experiment_info.txt"
echo "Fin: $(date)" >> "$RESULTS_BASE/experiment_info.txt"
echo "Tiempo total: ${hours}h ${minutes}m ${seconds}s" >> "$RESULTS_BASE/experiment_info.txt"
echo "Instancias completadas: $completed_instances" >> "$RESULTS_BASE/experiment_info.txt"

echo -e "\n${YELLOW}Tip: Usa el archivo summary.csv para analizar los resultados${NC}"