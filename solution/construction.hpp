#ifndef _FILO2_CONSTRUCTION_HPP_
#define _FILO2_CONSTRUCTION_HPP_

#include <vector>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <map>
#include <iostream>
#include <cstdlib> 
#include <cstring> 

#include "../instance/Instance.hpp"
#include "../solution/Solution.hpp"

extern "C" {
    // Desactivar warnings de Clang/GCC
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
    // #include "edgegen.h"  <-- ELIMINADO PARA EVITAR ERRORES DE LINKEO
    #include "macrorus.h"
    #include "linkern.h"
}

#define LK_RANDOM   (0)
#define LK_NEIGHBOR (1)
#define LK_GREEDY   (2)
#define LK_BORUVKA  (3)
#define LK_QBORUVKA (4)
#define KICK_TYPE   CC_LK_WALK_KICK

namespace cobra {

    struct PolarCoord {
        int customer_id;
        double rho;
        double phi;
        int demand;
    };

    inline double call_linkern_and_reorder(int n, int init, const Instance &instance, std::vector<int> &assignment) {
        if (n < 2) return 0.0; 
        if (n == 2) { 
             return 2.0 * instance.get_cost(assignment[init], instance.get_depot());
        }

        int n_tsp = n + 1; 
        
        CCdatagroup dat;
        CCrandstate rstate;
        
        CCutil_init_datagroup(&dat);
        
        int seed = (int)CCutil_real_zeit();
        CCutil_sprand(seed, &rstate);

        dat.x = (double *) CC_SAFE_MALLOC(n_tsp, double);
        dat.y = (double *) CC_SAFE_MALLOC(n_tsp, double);
        
        int depot_id = instance.get_depot();
        dat.x[0] = instance.get_x_coordinate(depot_id);
        dat.y[0] = instance.get_y_coordinate(depot_id);

        for (int i = 0; i < n; ++i) {
            int cust_id = assignment[init + i];
            dat.x[i+1] = instance.get_x_coordinate(cust_id);
            dat.y[i+1] = instance.get_y_coordinate(cust_id);
        }

        int norm = CC_EUCLIDEAN; 
        CCutil_dat_setnorm(&dat, norm);
        CCutil_dat_getnorm(&dat, &norm); 

        int *incycle = (int *) CC_SAFE_MALLOC(n_tsp, int);
        int *outcycle = (int *) CC_SAFE_MALLOC(n_tsp, int);
        double val = 0.0;
        bool success = false;

        // Variables para aristas
        int ecount = 0;
        int *elist = (int *) NULL;

        // --- PASO 1: Tour Inicial (Boruvka) ---
        CCkdtree localkt;
        // Nota: Boruvka internamente no usa aristas explícitas, usa geometría
        if (CCkdtree_build(&localkt, n_tsp, &dat, NULL, &rstate) == 0) {
             int tempcount;
             int *templist;
             CCkdtree_quadrant_k_nearest(&localkt, n_tsp, 2, &dat, NULL, 1, &tempcount, &templist, 1, &rstate);
             if (CCkdtree_qboruvka_tour(&localkt, n_tsp, &dat, incycle, &val, &rstate) == 0) {
                 success = true;
             }
             CC_IFFREE(templist, int);
             CCkdtree_free(&localkt);
        }

        // --- PASO 2: Generar Aristas (Grafo Completo / Clique) ---
        // [SOLUCIÓN DEFINITIVA]
        // En lugar de depender de edgegen, generamos todas las conexiones posibles.
        // Para rutas de camiones (ej: 20-50 nodos), esto es trivial y muy rápido.
        if (success) {
            ecount = (n_tsp * (n_tsp - 1)) / 2;
            elist = (int *) CC_SAFE_MALLOC(2 * ecount, int);
            
            int k = 0;
            for (int i = 0; i < n_tsp; ++i) {
                for (int j = i + 1; j < n_tsp; ++j) {
                    elist[2*k] = i;
                    elist[2*k+1] = j;
                    k++;
                }
            }
        }

        // --- PASO 3: Lin-Kernighan ---
        if (success) {
            int run_silently = 1;
            // Linkern ahora recibe un grafo completo, imposible que falle por falta de aristas
            if (CClinkern_tour(n_tsp, &dat, ecount, elist, 10000000, n_tsp, incycle, outcycle, &val, 
                            run_silently, -1.0, -1.0, NULL, KICK_TYPE, &rstate) != 0) {
                std::cerr << "Linkern failed internal error." << std::endl;
                success = false;
            }
        } else {
            for(int i=0; i<n_tsp; ++i) outcycle[i] = (i + 1) % n_tsp;
        }

        // --- PASO 4: Reordenar ---
        if (success) {
            std::vector<int> optimized_segment;
            optimized_segment.reserve(n);
            
            int depot_idx_in_tour = -1;
            for(int i=0; i<n_tsp; ++i) {
                if(outcycle[i] == 0) {
                    depot_idx_in_tour = i;
                    break;
                }
            }

            if (depot_idx_in_tour != -1) {
                for(int i=1; i<n_tsp; ++i) {
                    int idx = outcycle[(depot_idx_in_tour + i) % n_tsp];
                    if(idx > 0) {
                        optimized_segment.push_back(assignment[init + (idx - 1)]);
                    }
                }
                for(int i=0; i<n; ++i) {
                    assignment[init + i] = optimized_segment[i];
                }
            }
        }

        // Limpieza
        CC_IFFREE(incycle, int);
        CC_IFFREE(outcycle, int);
        CC_IFFREE(elist, int); 
        CCutil_freedatagroup(&dat);

        return val;
    }

    // ... (El resto de funciones: pack_in_trucks, compute_routes_and_build, sweep SÍGUEN IGUAL)
    inline void pack_in_trucks(const Instance &instance, 
                               const std::vector<PolarCoord> &toPack, 
                               std::vector<int> &assignment) {
        
        int remainingCapacity = instance.get_vehicle_capacity();
        int depotId = instance.get_depot(); 

        for (const auto &pc : toPack) {
            if (remainingCapacity < pc.demand) {
                assignment.push_back(depotId);
                remainingCapacity = instance.get_vehicle_capacity();
            }
            
            assert(remainingCapacity >= pc.demand);
            assignment.push_back(pc.customer_id);
            remainingCapacity -= pc.demand;
        }
        assignment.push_back(depotId);
    }

    inline void compute_routes_and_build(const Instance &instance, 
                                         std::vector<int> &assignment, 
                                         Solution &solution) {
        int depotId = instance.get_depot();
        int pointerToInit = 0;
        int n_customers_in_route = 0;
        
        for (unsigned int i = 0; i < assignment.size(); ++i) {
            if (assignment[i] == depotId) {
                if (n_customers_in_route > 0) {
                    call_linkern_and_reorder(n_customers_in_route, pointerToInit, instance, assignment);
                    
                    int first_cust = assignment[pointerToInit];
                    solution.build_one_customer_route</*record_action=*/false>(first_cust);
                    int current_route_idx = solution.get_route_index(first_cust);

                    for (int k = 1; k < n_customers_in_route; ++k) {
                        int next_cust = assignment[pointerToInit + k];
                        solution.insert_vertex_before</*record_action=*/false>(
                            current_route_idx, 
                            depotId, 
                            next_cust
                        );
                    }
                }
                pointerToInit = i + 1;
                n_customers_in_route = 0;
            } else {
                n_customers_in_route++;
            }
        }
    }

    inline void sweep(const Instance &instance, Solution &solution, double th) {
        solution.reset();
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
        
        std::sort(polar_coords.begin(), polar_coords.end(), 
                [](const PolarCoord& a, const PolarCoord& b) {
                    return a.phi < b.phi;
                });

        double threshold_distance = th * max_rho;
        std::vector<PolarCoord> inner, outer;
        for (const auto& pc : polar_coords) {
            if (pc.rho < threshold_distance) inner.push_back(pc);
            else outer.push_back(pc);
        }

        std::vector<int> assignment;
        if (!inner.empty()) pack_in_trucks(instance, inner, assignment);
        if (!outer.empty()) pack_in_trucks(instance, outer, assignment);

        compute_routes_and_build(instance, assignment, solution);
        assert(solution.is_feasible());
    }

}  // namespace cobra

#endif