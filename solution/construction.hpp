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

    struct PolarCoord {
        int customer_id;
        double rho;
        double phi;
        int demand;
    };

    // --- CONCORDE WRAPPER (Igual que antes) ---
    inline double call_linkern_and_reorder(int n, int init, const Instance &instance, std::vector<int> &assignment) {
        if (n < 2) return 0.0; 
        if (n == 2) return 2.0 * instance.get_cost(assignment[init], instance.get_depot());

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

        // Aristas explícitas (Grafo completo)
        int ecount = (n_tsp * (n_tsp - 1)) / 2;
        int *elist = (int *) CC_SAFE_MALLOC(2 * ecount, int);
        int k = 0;
        for (int i = 0; i < n_tsp; ++i) {
            for (int j = i + 1; j < n_tsp; ++j) {
                elist[2*k] = i; elist[2*k+1] = j; k++;
            }
        }

        // Boruvka + Linkern
        CCkdtree localkt;
        if (CCkdtree_build(&localkt, n_tsp, &dat, NULL, &rstate) == 0) {
             int tempcount, *templist;
             CCkdtree_quadrant_k_nearest(&localkt, n_tsp, 2, &dat, NULL, 1, &tempcount, &templist, 1, &rstate);
             CCkdtree_qboruvka_tour(&localkt, n_tsp, &dat, incycle, &val, &rstate);
             CC_IFFREE(templist, int);
             CCkdtree_free(&localkt);
        }

        int run_silently = 1;
        if (CClinkern_tour(n_tsp, &dat, ecount, elist, 10000000, n_tsp, incycle, outcycle, &val, 
                        run_silently, -1.0, -1.0, NULL, KICK_TYPE, &rstate) == 0) {
            success = true;
        } else {
            for(int i=0; i<n_tsp; ++i) outcycle[i] = (i + 1) % n_tsp; // Fallback
        }

        if (success) {
            std::vector<int> optimized_segment;
            optimized_segment.reserve(n);
            int depot_idx = -1;
            for(int i=0; i<n_tsp; ++i) if(outcycle[i]==0) { depot_idx=i; break; }

            if (depot_idx != -1) {
                for(int i=1; i<n_tsp; ++i) {
                    int idx = outcycle[(depot_idx + i) % n_tsp];
                    if(idx > 0) optimized_segment.push_back(assignment[init + (idx - 1)]);
                }
                for(int i=0; i<n; ++i) assignment[init + i] = optimized_segment[i];
            }
        }

        CC_IFFREE(incycle, int); CC_IFFREE(outcycle, int);
        CC_IFFREE(elist, int); CCutil_freedatagroup(&dat);
        return val;
    }

    // --- LÓGICA DE PADDING (NUEVO) ---
    // Ajusta el umbral para que el grupo 'outer' sea múltiplo de la capacidad
    inline double get_smart_threshold(const Instance &instance, 
                                    const std::vector<PolarCoord> &sorted_by_rho, 
                                    double initial_th_ratio) {
        
        int capacity = instance.get_vehicle_capacity();
        double max_rho = sorted_by_rho.back().rho;
        double target_rho = initial_th_ratio * max_rho;

        // Buscar el cliente más cercano a ese radio target
        auto it = std::lower_bound(sorted_by_rho.begin(), sorted_by_rho.end(), target_rho, 
            [](const PolarCoord &pc, double val) { return pc.rho < val; });

        // Calcular cuánta demanda hay desde 'it' hasta el final (grupo Outer)
        int outer_demand = 0;
        for (auto curr = it; curr != sorted_by_rho.end(); ++curr) {
            outer_demand += curr->demand;
        }

        // Ajustar 'it' para que outer_demand sea múltiplo de capacity
        // (Intentamos reducir el residuo)
        
        // Estrategia simple: Mover el corte hacia afuera hasta que encaje mejor
        // o hasta que se nos acabe el rango.
        auto best_it = it;
        int best_remainder = outer_demand % capacity;
        
        // Buscar en una ventana local (ej: mover el corte +/- 50 clientes)
        // para minimizar el desperdicio del último camión
        auto search_start = it;
        auto search_end = it;
        int range = 100; // mirar 100 clientes alrededor
        
        while(range > 0 && search_start != sorted_by_rho.begin()) { search_start--; range--; }
        range = 100;
        while(range > 0 && search_end != sorted_by_rho.end()) { search_end++; range--; }

        // Recalcular demanda para la ventana inicial
        int current_demand = 0;
        for (auto curr = search_start; curr != sorted_by_rho.end(); ++curr) {
            current_demand += curr->demand;
        }

        int min_waste = capacity; // Queremos minimizar (capacidad - residuo) o simplemente residuo bajo
        
        for (auto curr = search_start; curr != search_end; ++curr) {
            int remainder = current_demand % capacity;
            // El "desperdicio" es cuánto espacio sobra en el último camión.
            // Si remainder es 0, desperdicio 0. Si remainder es 1, desperdicio capacity - 1.
            // Queremos que el remainder sea cercano a 0 (camión lleno) o cercano a capacity (otro camión lleno).
            
            int waste = (remainder == 0) ? 0 : (capacity - remainder);

            if (waste < min_waste) {
                min_waste = waste;
                best_it = curr;
            }
            current_demand -= curr->demand; // Al avanzar el iterador, ese cliente pasa de Outer a Inner
        }

        if (best_it == sorted_by_rho.end()) return max_rho + 0.1;
        return best_it->rho;
    }

    inline void pack_in_trucks(const Instance &instance, const std::vector<PolarCoord> &toPack, std::vector<int> &assignment) {
        int remainingCapacity = instance.get_vehicle_capacity();
        int depotId = instance.get_depot(); 
        for (const auto &pc : toPack) {
            if (remainingCapacity < pc.demand) {
                assignment.push_back(depotId);
                remainingCapacity = instance.get_vehicle_capacity();
            }
            assignment.push_back(pc.customer_id);
            remainingCapacity -= pc.demand;
        }
        assignment.push_back(depotId);
    }

    inline void compute_routes_and_build(const Instance &instance, std::vector<int> &assignment, Solution &solution) {
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
                        solution.insert_vertex_before</*record_action=*/false>(current_route_idx, depotId, assignment[pointerToInit + k]);
                    }
                }
                pointerToInit = i + 1;
                n_customers_in_route = 0;
            } else {
                n_customers_in_route++;
            }
        }
    }

    // --- SWEEP PRINCIPAL ACTUALIZADO ---
    inline void sweep(const Instance &instance, Solution &solution, double th_ratio) {
        solution.reset();
        
        std::vector<PolarCoord> polar_coords;
        polar_coords.reserve(instance.get_customers_num());
        
        const auto depot = instance.get_depot();
        const auto depot_x = instance.get_x_coordinate(depot);
        const auto depot_y = instance.get_y_coordinate(depot);
        
        for (auto i = instance.get_customers_begin(); i < instance.get_customers_end(); i++) {
            const auto dx = instance.get_x_coordinate(i) - depot_x;
            const auto dy = instance.get_y_coordinate(i) - depot_y;
            PolarCoord pc;
            pc.customer_id = i;
            pc.rho = std::sqrt(dx * dx + dy * dy);
            pc.phi = std::atan2(dy, dx);
            pc.demand = instance.get_demand(i);
            polar_coords.push_back(pc);
        }
        
        // Copia para ordenar por Rho (para el Padding)
        std::vector<PolarCoord> sorted_by_rho = polar_coords;
        std::sort(sorted_by_rho.begin(), sorted_by_rho.end(), 
            [](const PolarCoord& a, const PolarCoord& b) { return a.rho < b.rho; });

        // [MEJORA] Calcular el umbral inteligente usando Padding
        double threshold_distance = get_smart_threshold(instance, sorted_by_rho, th_ratio);
        
        // Ordenar original por Phi (Ángulo) para el barrido
        std::sort(polar_coords.begin(), polar_coords.end(), 
                [](const PolarCoord& a, const PolarCoord& b) { return a.phi < b.phi; });

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

}

#endif