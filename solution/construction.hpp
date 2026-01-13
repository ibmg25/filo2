#ifndef _FILO2_CONSTRUCTION_HPP_
#define _FILO2_CONSTRUCTION_HPP_

#include <vector>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <cstring> // Para memset

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
        double rho;        // Distancia al depósito
        double phi;        // Ángulo respecto al eje X
        int demand;
    };
    
    inline double call_linkern_optimized(int n, int init, const Instance &instance, 
                                         std::vector<int> &assignment) {
        // ... (Inicio de la función igual que antes: verificaciones n<2, n=3, etc.) ...
        if (n < 2) return 0.0;
        int depot_id = instance.get_depot();
        
        if (n == 2) {
             return 2.0 * (instance.get_cost(depot_id, assignment[init]) + 
                         instance.get_cost(assignment[init], assignment[init + 1]));
        }
        
        if (n == 3) {
             int c1 = assignment[init];
             int c2 = assignment[init + 1];
             int c3 = assignment[init + 2];
             return instance.get_cost(depot_id, c1) + 
                    instance.get_cost(c1, c2) + 
                    instance.get_cost(c2, c3) + 
                    instance.get_cost(c3, depot_id);
        }

        int n_tsp = n + 1;
        
        // ... (Buffers estáticos igual que antes) ...
        static thread_local std::vector<double> x_buf;
        static thread_local std::vector<double> y_buf;
        static thread_local std::vector<int> incycle_buf;
        static thread_local std::vector<int> outcycle_buf;
        static thread_local std::vector<int> elist_buf;

        if (x_buf.size() < (size_t)n_tsp) {
            x_buf.resize(n_tsp);
            y_buf.resize(n_tsp);
            incycle_buf.resize(n_tsp);
            outcycle_buf.resize(n_tsp);
        }
        
        // --- CAMBIO 1: AUMENTAR VENTANA ---
        // Subimos a 40 para capturar más geometría local.
        int neighbors_window = 40; 
        size_t estimated_edges = n_tsp * (neighbors_window + 5); 
        if (elist_buf.size() < 2 * estimated_edges) {
            elist_buf.resize(2 * estimated_edges);
        }

        // ... (Setup de CCdatagroup y coordenadas igual que antes) ...
        CCdatagroup dat;
        std::memset(&dat, 0, sizeof(CCdatagroup)); 
        dat.x = x_buf.data();
        dat.y = y_buf.data();
        dat.x[0] = instance.get_x_coordinate(depot_id);
        dat.y[0] = instance.get_y_coordinate(depot_id);
        for (int i = 0; i < n; ++i) {
            int cust_id = assignment[init + i];
            dat.x[i + 1] = instance.get_x_coordinate(cust_id);
            dat.y[i + 1] = instance.get_y_coordinate(cust_id);
        }
        CCutil_dat_setnorm(&dat, CC_EUCLIDEAN);
        
        CCrandstate rstate;
        int seed = 1234 + init; 
        CCutil_sprand(seed, &rstate);

        int* incycle = incycle_buf.data();
        for(int i=0; i<n_tsp; ++i) incycle[i] = i;

        int* outcycle = outcycle_buf.data();
        double val = 0.0;
        
        // Generación de aristas (Mismo código, solo cambia el efecto de neighbors_window)
        int* elist = elist_buf.data();
        int k_idx = 0;
        for (int i = 0; i < n_tsp; ++i) {
            for (int offset = 1; offset <= neighbors_window; ++offset) {
                int j = i + offset;
                if (j < n_tsp) {
                    elist[2 * k_idx] = i;
                    elist[2 * k_idx + 1] = j;
                    k_idx++;
                }
            }
            if (i > 0 && i > neighbors_window) { 
                elist[2 * k_idx] = 0;
                elist[2 * k_idx + 1] = i;
                k_idx++;
            }
        }
        int ecount = k_idx;

        // --- CAMBIO 2: REACTIVAR KICKS ---
        int run_silently = 1;
        // Damos un poco de "juego" para escapar de óptimos locales malos
        int max_kicks = 20; 

        bool lk_success = (CClinkern_tour(n_tsp, &dat, ecount, elist, max_kicks, 
                                        0, incycle, outcycle, &val, 
                                        run_silently, -1.0, -1.0, NULL, KICK_TYPE, &rstate) == 0);

        // ... (Reconstrucción igual que antes) ...
        if (lk_success) {
            int depot_idx = -1;
            for (int i = 0; i < n_tsp; ++i) {
                if (outcycle[i] == 0) { depot_idx = i; break; }
            }
            if (depot_idx != -1) {
                static thread_local std::vector<int> opt_segment;
                if (opt_segment.size() < (size_t)n) opt_segment.resize(n);
                for (int i = 1; i < n_tsp; ++i) {
                    int idx = outcycle[(depot_idx + i) % n_tsp];
                    if (idx > 0) opt_segment[i-1] = assignment[init + (idx - 1)];
                }
                for (int i = 0; i < n; ++i) assignment[init + i] = opt_segment[i];
            }
        }
        
        return val;
    }
    
    inline void pack_in_trucks(const Instance &instance, 
                               const std::vector<PolarCoord> &toPack, 
                               std::vector<int> &assignment) {
        
        int remaining_capacity = instance.get_vehicle_capacity();
        int depot_id = instance.get_depot();
        
        for (const auto &pc : toPack) {
            if (remaining_capacity < pc.demand) {
                assignment.push_back(depot_id); 
                remaining_capacity = instance.get_vehicle_capacity();
            }
            assignment.push_back(pc.customer_id);
            remaining_capacity -= pc.demand;
            assert(remaining_capacity >= 0);
        }
        assignment.push_back(depot_id);
    }
    
    inline void compute_routes_and_build(const Instance &instance, 
                                         std::vector<int> &assignment, 
                                         Solution &solution) {
        
        int depot_id = instance.get_depot();
        int pointer_to_init = 0;
        int n_customers_in_route = 0;
        
        for (unsigned int i = 0; i < assignment.size(); ++i) {
            if (assignment[i] == depot_id) {
                if (n_customers_in_route > 0) {
                    call_linkern_optimized(n_customers_in_route, pointer_to_init, 
                                          instance, assignment);
                    
                    int first_cust = assignment[pointer_to_init];
                    solution.build_one_customer_route</*record_action=*/false>(first_cust);
                    int current_route_idx = solution.get_route_index(first_cust);
                    
                    for (int k = 1; k < n_customers_in_route; ++k) {
                        solution.insert_vertex_before</*record_action=*/false>(
                            current_route_idx, depot_id, assignment[pointer_to_init + k]);
                    }
                }
                pointer_to_init = i + 1;
                n_customers_in_route = 0;
            } else {
                n_customers_in_route++;
            }
        }
    }
    
    // --- SWEEP MODIFICADO: AHORA ACEPTA ROTATION_IDX ---
    inline void sweep(const Instance &instance, Solution &solution, double th_ratio, int rotation_idx = 0) {
        
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
        
        double threshold_distance = th_ratio * max_rho;
        
        std::sort(polar_coords.begin(), polar_coords.end(), 
                  [](const PolarCoord& a, const PolarCoord& b) { 
                      return a.phi < b.phi; 
                  });
        
        // --- AQUI SE APLICA LA ROTACIÓN ---
        if (rotation_idx > 0 && rotation_idx < (int)polar_coords.size()) {
            std::rotate(polar_coords.begin(), polar_coords.begin() + rotation_idx, polar_coords.end());
        }
        
        std::vector<PolarCoord> inner, outer;
        for (const auto& pc : polar_coords) {
            if (pc.rho < threshold_distance) {
                inner.push_back(pc);
            } else {
                outer.push_back(pc);
            }
        }
        
        std::vector<int> assignment;
        if (!inner.empty()) pack_in_trucks(instance, inner, assignment);
        if (!outer.empty()) pack_in_trucks(instance, outer, assignment);
        
        compute_routes_and_build(instance, assignment, solution);
        assert(solution.is_feasible());
    }

} // namespace cobra

#endif