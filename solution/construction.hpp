#ifndef _FILO2_CONSTRUCTION_HPP_
#define _FILO2_CONSTRUCTION_HPP_

#include <vector>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <iostream>

#include "../instance/Instance.hpp"
#include "../solution/Solution.hpp"

extern "C" {
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wkeyword-macro"
    #define new _new_var_
    #define class _class_var_
    #ifdef NDEBUG
    #undef NDEBUG
    #endif
    #include "config.h"
    #include "concorde.h"
    #undef new
    #undef class
    #pragma clang diagnostic pop 
    #include "machdefs.h"
    #include "util.h"
    #include "kdtree.h"
    #include "macrorus.h"
    #include "linkern.h"
}

#define KICK_TYPE CC_LK_WALK_KICK

namespace cobra {

    // =========================================================================
    // ESTRUCTURA DE COORDENADAS POLARES
    // =========================================================================
    
    struct PolarCoord {
        int customer_id;
        double rho;        // Distancia al depósito
        double phi;        // Ángulo respecto al eje X
        int demand;
    };

    // =========================================================================
    // WRAPPER OPTIMIZADO DE CONCORDE LIN-KERNIGHAN
    // =========================================================================
    
    inline double call_linkern_optimized(int n, int init, const Instance &instance, 
                                         std::vector<int> &assignment) {
        // Casos triviales - evitar overhead de Concorde
        if (n < 2) return 0.0;
        
        int depot_id = instance.get_depot();
        
        if (n == 2) {
            // Ruta: depot -> c1 -> c2 -> depot
            return 2.0 * (instance.get_cost(depot_id, assignment[init]) + 
                         instance.get_cost(assignment[init], assignment[init + 1]));
        }
        
        if (n == 3) {
            // Ruta: depot -> c1 -> c2 -> c3 -> depot
            int c1 = assignment[init];
            int c2 = assignment[init + 1];
            int c3 = assignment[init + 2];
            return instance.get_cost(depot_id, c1) + 
                   instance.get_cost(c1, c2) + 
                   instance.get_cost(c2, c3) + 
                   instance.get_cost(c3, depot_id);
        }

        // =====================================================================
        // Preparar datos para Concorde (incluye depósito)
        // =====================================================================
        
        int n_tsp = n + 1; // n clientes + depósito
        CCdatagroup dat;
        CCrandstate rstate;
        CCutil_init_datagroup(&dat);
        
        int seed = (int)CCutil_real_zeit();
        CCutil_sprand(seed, &rstate);

        dat.x = (double *) CC_SAFE_MALLOC(n_tsp, double);
        dat.y = (double *) CC_SAFE_MALLOC(n_tsp, double);
        
        if (!dat.x || !dat.y) {
            std::cerr << "ERROR: Memory allocation failed in Concorde wrapper\n";
            CCutil_freedatagroup(&dat);
            return INFINITY;
        }
        
        // Índice 0: depósito
        dat.x[0] = instance.get_x_coordinate(depot_id);
        dat.y[0] = instance.get_y_coordinate(depot_id);
        
        // Índices 1..n: clientes
        for (int i = 0; i < n; ++i) {
            int cust_id = assignment[init + i];
            dat.x[i + 1] = instance.get_x_coordinate(cust_id);
            dat.y[i + 1] = instance.get_y_coordinate(cust_id);
        }

        int norm = CC_EUCLIDEAN;
        CCutil_dat_setnorm(&dat, norm);

        // =====================================================================
        // Crear solución inicial con Boruvka + K-d tree
        // =====================================================================
        
        int *incycle = (int *) CC_SAFE_MALLOC(n_tsp, int);
        int *outcycle = (int *) CC_SAFE_MALLOC(n_tsp, int);
        double val = 0.0;

        if (!incycle || !outcycle) {
            std::cerr << "ERROR: Memory allocation failed for cycles\n";
            CC_IFFREE(incycle, int);
            CC_IFFREE(outcycle, int);
            CCutil_freedatagroup(&dat);
            return INFINITY;
        }

        CCkdtree localkt;
        bool kdtree_built = false;
        bool boruvka_ok = false;
        
        // Construir K-d tree
        if (CCkdtree_build(&localkt, n_tsp, &dat, NULL, &rstate) == 0) {
            kdtree_built = true;
            
            // Generar vecinos cercanos para heurística inicial
            int tempcount, *templist = NULL;
            if (CCkdtree_quadrant_k_nearest(&localkt, n_tsp, 2, &dat, NULL, 1, 
                                           &tempcount, &templist, 1, &rstate) == 0) {
                
                // Solución inicial con Boruvka
                if (CCkdtree_qboruvka_tour(&localkt, n_tsp, &dat, incycle, &val, &rstate) == 0) {
                    boruvka_ok = true;
                }
                
                CC_IFFREE(templist, int);
            }
        }

        // Fallback: tour trivial si Boruvka falla
        if (!boruvka_ok) {
            for (int i = 0; i < n_tsp; ++i) {
                incycle[i] = i;
            }
        }

        // =====================================================================
        // Aplicar Lin-Kernighan para optimización
        // =====================================================================
        
        // Generar lista de arcos para Lin-Kernighan (grafo completo)
        int ecount = (n_tsp * (n_tsp - 1)) / 2;
        int *elist = (int *) CC_SAFE_MALLOC(2 * ecount, int);
        
        if (!elist) {
            std::cerr << "ERROR: Memory allocation failed for edge list\n";
            CC_IFFREE(incycle, int);
            CC_IFFREE(outcycle, int);
            if (kdtree_built) CCkdtree_free(&localkt);
            CCutil_freedatagroup(&dat);
            return INFINITY;
        }
        
        int k = 0;
        for (int i = 0; i < n_tsp; ++i) {
            for (int j = i + 1; j < n_tsp; ++j) {
                elist[2 * k] = i;
                elist[2 * k + 1] = j;
                k++;
            }
        }

        int run_silently = 1;
        int in_repeater = std::min(n_tsp, 5); // Parámetro de intensidad LK
        
        bool lk_success = false;
        if (CClinkern_tour(n_tsp, &dat, ecount, elist, 10000000, in_repeater, 
                          incycle, outcycle, &val, run_silently, -1.0, -1.0, 
                          NULL, KICK_TYPE, &rstate) == 0) {
            lk_success = true;
        }

        // =====================================================================
        // Reordenar assignment según tour óptimo encontrado
        // =====================================================================
        
        if (lk_success) {
            // Encontrar posición del depósito en el tour
            int depot_idx = -1;
            for (int i = 0; i < n_tsp; ++i) {
                if (outcycle[i] == 0) { // 0 = índice del depósito
                    depot_idx = i;
                    break;
                }
            }

            if (depot_idx != -1) {
                std::vector<int> optimized_segment;
                optimized_segment.reserve(n);
                
                // Recorrer tour desde depósito, agregando solo clientes
                for (int i = 1; i < n_tsp; ++i) {
                    int idx = outcycle[(depot_idx + i) % n_tsp];
                    if (idx > 0) { // Ignorar depósito (idx=0)
                        optimized_segment.push_back(assignment[init + (idx - 1)]);
                    }
                }
                
                // Actualizar assignment con orden optimizado
                for (int i = 0; i < n; ++i) {
                    assignment[init + i] = optimized_segment[i];
                }
            }
        } else {
            // Si Lin-Kernighan falló, mantener orden original
            std::cerr << "WARNING: Lin-Kernighan failed, keeping original order\n";
        }

        // =====================================================================
        // Limpieza de memoria
        // =====================================================================
        
        CC_IFFREE(incycle, int);
        CC_IFFREE(outcycle, int);
        CC_IFFREE(elist, int);
        if (kdtree_built) CCkdtree_free(&localkt);
        CCutil_freedatagroup(&dat);

        return val;
    }

    // =========================================================================
    // EMPAQUETADO DE CLIENTES EN CAMIONES (ORDEN ANGULAR ESTRICTO)
    // =========================================================================
    
    inline void pack_in_trucks(const Instance &instance, 
                               const std::vector<PolarCoord> &toPack, 
                               std::vector<int> &assignment) {
        
        int remaining_capacity = instance.get_vehicle_capacity();
        int depot_id = instance.get_depot();
        
        for (const auto &pc : toPack) {
            // Si el cliente no cabe en el camión actual, cerrar ruta y abrir nueva
            if (remaining_capacity < pc.demand) {
                assignment.push_back(depot_id); // Marca de fin de ruta
                remaining_capacity = instance.get_vehicle_capacity();
            }
            
            // Agregar cliente a la ruta actual
            assignment.push_back(pc.customer_id);
            remaining_capacity -= pc.demand;
            
            assert(remaining_capacity >= 0);
        }
        
        // Cerrar última ruta
        assignment.push_back(depot_id);
    }

    // =========================================================================
    // CONSTRUCCIÓN DE RUTAS CON CONCORDE + INSERCIÓN EN SOLUCIÓN
    // =========================================================================
    
    inline void compute_routes_and_build(const Instance &instance, 
                                         std::vector<int> &assignment, 
                                         Solution &solution) {
        
        int depot_id = instance.get_depot();
        int pointer_to_init = 0;
        int n_customers_in_route = 0;
        
        for (unsigned int i = 0; i < assignment.size(); ++i) {
            if (assignment[i] == depot_id) {
                // Fin de ruta detectado
                
                if (n_customers_in_route > 0) {
                    // Optimizar ruta con Concorde
                    call_linkern_optimized(n_customers_in_route, pointer_to_init, 
                                          instance, assignment);
                    
                    // Construir ruta en la solución
                    int first_cust = assignment[pointer_to_init];
                    solution.build_one_customer_route</*record_action=*/false>(first_cust);
                    int current_route_idx = solution.get_route_index(first_cust);
                    
                    // Insertar resto de clientes en orden
                    for (int k = 1; k < n_customers_in_route; ++k) {
                        solution.insert_vertex_before</*record_action=*/false>(
                            current_route_idx, depot_id, assignment[pointer_to_init + k]);
                    }
                }
                
                // Preparar para siguiente ruta
                pointer_to_init = i + 1;
                n_customers_in_route = 0;
                
            } else {
                // Cliente en ruta actual
                n_customers_in_route++;
            }
        }
    }

    // =========================================================================
    // ALGORITMO SWEEP PRINCIPAL
    // =========================================================================
    
    inline void sweep(const Instance &instance, Solution &solution, double th_ratio) {
        
        solution.reset();
        
        // =====================================================================
        // 1. Calcular coordenadas polares de todos los clientes
        // =====================================================================
        
        std::vector<PolarCoord> polar_coords;
        polar_coords.reserve(instance.get_customers_num());
        
        const auto depot = instance.get_depot();
        const auto depot_x = instance.get_x_coordinate(depot);
        const auto depot_y = instance.get_y_coordinate(depot);
        
        double max_rho = 0.0;
        
        for (auto i = instance.get_customers_begin(); i < instance.get_customers_end(); i++) {
            const auto dx = instance.get_x_coordinate(i) - depot_x;
            const auto dy = instance.get_y_coordinate(i) - depot_y;
            
            PolarCoord pc;
            pc.customer_id = i;
            pc.rho = std::sqrt(dx * dx + dy * dy);
            pc.phi = std::atan2(dy, dx);
            pc.demand = instance.get_demand(i);
            
            polar_coords.push_back(pc);
            
            max_rho = std::max(max_rho, pc.rho);
        }
        
        // =====================================================================
        // 2. Calcular threshold de distancia radial
        // =====================================================================
        
        double threshold_distance = th_ratio * max_rho;
        
        // =====================================================================
        // 3. Ordenar clientes por ángulo (phi) - ORDEN ANGULAR SWEEP
        // =====================================================================
        
        std::sort(polar_coords.begin(), polar_coords.end(), 
                  [](const PolarCoord& a, const PolarCoord& b) { 
                      return a.phi < b.phi; 
                  });

        // =====================================================================
        // 4. Particionar en INNER (cerca) y OUTER (lejos)
        // =====================================================================
        
        std::vector<PolarCoord> inner, outer;
        
        for (const auto& pc : polar_coords) {
            if (pc.rho < threshold_distance) {
                inner.push_back(pc);
            } else {
                outer.push_back(pc);
            }
        }

        // =====================================================================
        // 5. Empaquetar clientes en camiones (orden angular estricto)
        // =====================================================================
        
        std::vector<int> assignment;
        
        if (!inner.empty()) {
            pack_in_trucks(instance, inner, assignment);
        }
        
        if (!outer.empty()) {
            pack_in_trucks(instance, outer, assignment);
        }

        // =====================================================================
        // 6. Optimizar cada ruta con Concorde y construir solución final
        // =====================================================================
        
        compute_routes_and_build(instance, assignment, solution);
        
        // Verificación final
        assert(solution.is_feasible());
    }

} // namespace cobra

#endif