#include "construction.hpp"
#include "../base/Timer.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <iostream>
#include <cstring> 

// --- FIX: Compatibilidad C++ con Concorde ---
#define new new_
#define class class_
extern "C" {
    #include "concorde.h"
}
#undef new
#undef class
// ---------------------------------------------

namespace cobra {

    struct PolarNode {
        int id;
        double angle;
        double rho;
    };

    double optimize_route_with_linkern(const Instance &instance, std::vector<int>& route_nodes, CCrandstate& rstate) {
        
        // Caso trivial
        if (route_nodes.size() <= 1) {
            double cost = 0.0;
            int prev = instance.get_depot();
            for (int node : route_nodes) {
                cost += instance.get_cost(prev, node);
                prev = node;
            }
            cost += instance.get_cost(prev, instance.get_depot());
            return cost;
        }

        int n = route_nodes.size(); 
        int n_total = n + 1;        

        CCdatagroup dat;
        CCutil_init_datagroup(&dat);

        dat.x = (double *) malloc(n_total * sizeof(double));
        dat.y = (double *) malloc(n_total * sizeof(double));

        if (!dat.x || !dat.y) {
            std::cerr << "[FILO2] Out of memory for coords\n";
            if (dat.x) free(dat.x);
            if (dat.y) free(dat.y);
            return 1e9;
        }

        std::vector<int> map_idx(n_total);
        
        // Nodo 0 = Depósito
        dat.x[0] = instance.get_x_coordinate(instance.get_depot());
        dat.y[0] = instance.get_y_coordinate(instance.get_depot());
        map_idx[0] = instance.get_depot();

        for(int i = 0; i < n; ++i) {
            int customer = route_nodes[i];
            dat.x[i+1] = instance.get_x_coordinate(customer);
            dat.y[i+1] = instance.get_y_coordinate(customer);
            map_idx[i+1] = customer;
        }

        CCutil_dat_setnorm(&dat, CC_EUCLIDEAN); 

        // outcycle contendrá la permutación de nodos en el tour
        std::vector<int> outcycle(n_total, -1);
        double val = 0.0;
        
        // Linkern parameters
        int run_silent = 1; 

        int res = CClinkern_tour(
            n_total,                
            &dat,                   
            0, NULL,                
            10000000,               
            0,                      
            NULL,                   
            outcycle.data(),        
            &val,                   
            run_silent,             
            0.0, 0.0,               
            NULL,                   
            CC_LK_WALK_KICK,        
            &rstate                 
        );

        if (res != 0) {
            CCutil_freedatagroup(&dat);
            return 1e9;
        }

        // --- CORRECCIÓN CRÍTICA: Interpretación de outcycle como Permutación ---
        
        // 1. Buscar la posición del depósito (nodo 0 de Concorde) en el tour
        int depot_pos = -1;
        for (int i = 0; i < n_total; ++i) {
            if (outcycle[i] == 0) {
                depot_pos = i;
                break;
            }
        }

        if (depot_pos == -1) {
            // El depósito no está en el tour (algo muy raro)
            CCutil_freedatagroup(&dat);
            return 1e9;
        }

        // 2. Reconstruir la ruta rotando para empezar después del depósito
        std::vector<int> optimized_nodes;
        optimized_nodes.reserve(n);
        
        double cost = 0.0;
        int prev_filo_node = instance.get_depot();

        for(int k = 1; k < n_total; ++k) {
            // Recorrer circularmente: (depot_pos + k) % n_total
            int idx_in_outcycle = (depot_pos + k) % n_total;
            int concorde_node = outcycle[idx_in_outcycle];
            
            // Mapear de vuelta a ID de FILO2
            int filo_node = map_idx[concorde_node];
            
            // Sanity check
            if (filo_node == instance.get_depot()) {
                 // Esto no debería pasar si la lógica del bucle es k=1..n
                 continue; 
            }

            optimized_nodes.push_back(filo_node);
            cost += instance.get_cost(prev_filo_node, filo_node);
            prev_filo_node = filo_node;
        }
        
        // Cerrar el ciclo
        cost += instance.get_cost(prev_filo_node, instance.get_depot());
        
        // Actualizar la ruta original
        route_nodes = std::move(optimized_nodes);

        CCutil_freedatagroup(&dat); 
        return cost;
    }

    void sweep_concorde_heuristic(const Instance &instance, Solution &solution, int seed) {
        
        solution.reset();
        
        int num_customers = instance.get_customers_num();
        std::vector<PolarNode> polar_nodes;
        polar_nodes.reserve(num_customers);

        double depot_x = instance.get_x_coordinate(instance.get_depot());
        double depot_y = instance.get_y_coordinate(instance.get_depot());
        double max_rho = 0.0;

        for (int i = instance.get_customers_begin(); i < instance.get_customers_end(); ++i) {
            double dx = instance.get_x_coordinate(i) - depot_x;
            double dy = instance.get_y_coordinate(i) - depot_y;
            double rho = std::hypot(dx, dy);
            double angle = std::atan2(dy, dx);
            
            if (rho > max_rho) max_rho = rho;
            polar_nodes.push_back({i, angle, rho});
        }

        std::sort(polar_nodes.begin(), polar_nodes.end(), [](const PolarNode& a, const PolarNode& b) {
            return a.angle < b.angle;
        });

        CCrandstate rstate;
        CCutil_sprand(seed, &rstate);

        double best_global_cost = std::numeric_limits<double>::max();
        std::vector<std::vector<int>> best_routes_config;
        bool found_solution = false;

        // Usamos 100 pasos para barrer bien el threshold
        int steps = 100; 
        
        for (int t = 0; t < steps; ++t) {
            double threshold = t / (double)steps;
            double rho_limit = max_rho * threshold;

            std::vector<int> visit_order;
            visit_order.reserve(num_customers);

            for (const auto& node : polar_nodes) {
                if (node.rho <= rho_limit) visit_order.push_back(node.id);
            }
            for (const auto& node : polar_nodes) {
                if (node.rho > rho_limit) visit_order.push_back(node.id);
            }

            std::vector<std::vector<int>> current_routes;
            std::vector<int> current_route;
            int current_load = 0;
            double current_total_cost = 0.0;
            int capacity = instance.get_vehicle_capacity();
            bool failure = false;

            auto process_route = [&](std::vector<int>& route) {
                if (route.empty()) return;
                double r_cost = optimize_route_with_linkern(instance, route, rstate);
                if (r_cost >= 1e8) { 
                    failure = true;
                }
                current_total_cost += r_cost;
                current_routes.push_back(route);
            };

            for (int customer : visit_order) {
                int demand = instance.get_demand(customer);
                if (current_load + demand > capacity) {
                    process_route(current_route);
                    if (failure) break;
                    current_route.clear();
                    current_route.push_back(customer);
                    current_load = demand;
                } else {
                    current_route.push_back(customer);
                    current_load += demand;
                }
            }
            if (!failure) process_route(current_route);

            if (!failure && current_total_cost < best_global_cost) {
                best_global_cost = current_total_cost;
                best_routes_config = std::move(current_routes);
                found_solution = true;
            }
        }

        if (!found_solution) {
            std::cerr << "[FILO2] Error: Sweep heuristic failed to find any valid solution.\n";
            // Fallback trivial
            for (int i = instance.get_customers_begin(); i < instance.get_customers_end(); ++i) {
                solution.build_one_customer_route(i);
            }
            return;
        }

        // Construcción en FILO2
        for (const auto& route_vec : best_routes_config) {
            if (route_vec.empty()) continue;
            int route_idx = solution.build_one_customer_route(route_vec[0]);
            for (size_t i = 1; i < route_vec.size(); ++i) {
                solution.insert_vertex_before(route_idx, instance.get_depot(), route_vec[i]);
            }
        }
        
        if (!solution.is_feasible(true, true)) { 
             std::cerr << "[FILO2] CRITICAL ERROR: Constructed solution is infeasible!\n";
             exit(1);
        } else {
            #ifdef VERBOSE
            std::cout << "Sweep constructed feasible solution with cost: " << solution.get_cost() << "\n";
            #endif
        }
    }

}