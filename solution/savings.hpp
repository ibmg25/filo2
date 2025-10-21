#ifndef _FILO2_SOLUTIONALGORITHMS_HPP_
#define _FILO2_SOLUTIONALGORITHMS_HPP_

#include "../base/Timer.hpp"
#include "Solution.hpp"


namespace cobra {

    // Limited savings algorithm.
    inline void clarke_and_wright(const Instance &instance, Solution &solution, const double lambda, int neighbors_num) {

        solution.reset();

        for (auto i = instance.get_customers_begin(); i < instance.get_customers_end(); i++) {
            solution.build_one_customer_route</*record_acion=*/false>(i);
        }
        assert(solution.is_feasible());

        neighbors_num = std::min(instance.get_customers_num() - 1, neighbors_num);

        const auto savings_num = instance.get_customers_num() * neighbors_num;

        struct Saving {
            int i;
            int j;
            double value;
        };

        auto savings = std::vector<Saving>();
        savings.reserve(static_cast<unsigned long>(savings_num));

        for (auto i = instance.get_customers_begin(); i < instance.get_customers_end(); i++) {

            for (auto n = 1u, added = 0u; added < static_cast<unsigned int>(neighbors_num) && n < instance.get_neighbors_of(i).size();
                 n++) {

                const auto j = instance.get_neighbors_of(i)[n];

                if (i < j) {

                    const double value = +instance.get_cost(i, instance.get_depot()) + instance.get_cost(instance.get_depot(), j) -
                                         lambda * instance.get_cost(i, j);

                    savings.push_back({i, j, value});

                    added++;
                }
            }
        }


        std::sort(savings.begin(), savings.end(), [](const Saving &a, const Saving &b) { return a.value > b.value; });

#ifdef VERBOSE
        Timer timer;
#endif

        for (auto n = 0; n < static_cast<int>(savings.size()); ++n) {

            const auto &saving = savings[n];

            const auto i = saving.i;
            const auto j = saving.j;

            const auto iRoute = solution.get_route_index(i);
            const auto jRoute = solution.get_route_index(j);

            if (iRoute == jRoute) {
                continue;
            }

            if (solution.get_last_customer(iRoute) == i && solution.get_first_customer(jRoute) == j &&
                solution.get_route_load(iRoute) + solution.get_route_load(jRoute) <= instance.get_vehicle_capacity()) {

                solution.append_route(iRoute, jRoute);


            } else if (solution.get_last_customer(jRoute) == j && solution.get_first_customer(iRoute) == i &&
                       solution.get_route_load(iRoute) + solution.get_route_load(jRoute) <= instance.get_vehicle_capacity()) {

                solution.append_route(jRoute, iRoute);
            }

#ifdef VERBOSE
            if (timer.elapsed_time<std::chrono::seconds>() > 2) {
                std::cout << "Progress: " << 100.0 * (n + 1) / savings.size() << "%, Solution cost: " << solution.get_cost() << " \n";
                timer.reset();
            }
#endif
        }
        assert(solution.is_feasible());
    }

    inline void sweep(const Instance &instance, Solution &solution,
                    const double radius_threshold = 0.5,
                    const double radius_min = 0.5,
                    const double radius_max = 0.5,
                    const unsigned seed = 0) {
        solution.reset();
        
        struct PolarCoord {
            int customer_id;
            double rho;
            double phi;
            int demand;
        };
        
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

        // Preparar RNG: si seed == 0, usar random_device (no determinista), si seed != 0, determinista
        unsigned use_seed = seed;
        if (use_seed == 0) {
            std::random_device rd;
            use_seed = rd();
        }
        std::mt19937 rng(use_seed);

        // Determinar factor radial (radius_factor en [0,1])
        double radius_factor = radius_threshold; // default: usa el radius_threshold simple
        // Si el usuario pasó un rango válido y distinto, tomar aleatorio en [radius_min, radius_max]
        if (radius_min != radius_max) {
            double rmin = std::min(std::max(radius_min, 0.0), 1.0);
            double rmax = std::min(std::max(radius_max, 0.0), 1.0);
            if (rmin > rmax) std::swap(rmin, rmax);
            std::uniform_real_distribution<double> dist_r(rmin, rmax);
            radius_factor = dist_r(rng);
        } else {
            // si no hay rango, usar el radius_threshold (clamped)
            radius_factor = std::min(std::max(radius_threshold, 0.0), 1.0);
        }

        const double threshold_distance = radius_factor * max_rho;

        std::vector<PolarCoord> inner, outer;
        for (const auto& pc : polar_coords) {
            if (pc.rho < threshold_distance) {
                inner.push_back(pc);
            } else {
                outer.push_back(pc);
            }
        }

        auto assign_group = [&](const std::vector<PolarCoord>& group) {
            if (group.empty()) return;
            
            // elegir un cliente de inicio aleatorio dentro del grupo
            std::uniform_int_distribution<size_t> dist_idx(0, group.size() - 1);
            const size_t start_idx = dist_idx(rng);
            
            // crear la primera ruta con el cliente seleccionado
            const auto &first_pc = group[start_idx];
            solution.build_one_customer_route</*record_action=*/false>(first_pc.customer_id);
            int current_route_idx = solution.get_route_index(first_pc.customer_id);
            int current_load = first_pc.demand;
            
            // recorrer el resto circularmente
            for (size_t offset = 1; offset < group.size(); offset++) {
                const size_t idx = (start_idx + offset) % group.size();
                const auto &pc = group[idx];
                
                if (current_load + pc.demand > instance.get_vehicle_capacity()) {
                    // No cabe, crear nueva ruta
                    solution.build_one_customer_route</*record_action=*/false>(pc.customer_id);
                    current_route_idx = solution.get_route_index(pc.customer_id);
                    current_load = pc.demand;
                } else {
                    // Cabe, fusionar con ruta actual
                    solution.build_one_customer_route</*record_action=*/false>(pc.customer_id);
                    const int temp_route_idx = solution.get_route_index(pc.customer_id);
                    solution.append_route(current_route_idx, temp_route_idx);
                    current_load += pc.demand;
                }
            }
        };

        assign_group(inner);
        assign_group(outer);
        
        assert(solution.is_feasible());
    }

}  // namespace cobra

#endif