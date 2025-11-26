#ifndef _FILO2_INSTANCE_HPP_
#define _FILO2_INSTANCE_HPP_

#include <iostream>
#include <cassert>
#include <cmath>
#include <optional>
#include <string>

#include "../base/NonCopyable.hpp"
#include "Parser.hpp"

namespace cobra {
    namespace {
        inline double fastround(double value) {
            return static_cast<int>(value + 0.5);
        }
    }  // namespace

    // Manages a CVRP instance by providing a set of methods to query its properties.
    class Instance : private NonCopyable<Instance> {
    public:
        // Returns an optional containing a properly built instance if the parsing of the input file completes correctly, nullopt otherwise.
        // The parameter `num_neighbors` specifies the number of neighbors that are precomputed for each vertice. These neighbors are then
        // accessibly by using the `get_neighbors_of` method.
        static std::optional<Instance> make(const std::string& filepath, int num_neighbors);

        // Returns the instance size.
        inline int get_vertices_num() const {
            return demands.size();
        }

        // Returns the depot's index.
        inline int get_depot() const {
            return 0;
        }

        // Returns the vehicle capacity.
        inline int get_vehicle_capacity() const {
            return vehicle_capacity;
        };

        // Returns the number of customers.
        inline int get_customers_num() const {
            return get_vertices_num() - 1;
        };

        // Returns the index of the first customer.
        inline int get_customers_begin() const {
            return 1;
        }

        // Returns the index after the last customer.
        inline int get_customers_end() const {
            return get_vertices_num();
        }

        // Returns the index of the first vertex.
        inline int get_vertices_begin() const {
            return get_depot();
        }

        // Returns the index after the last vertex.
        inline int get_vertices_end() const {
            return get_customers_end();
        }

        // Returns the cost of arc (i, j).
        inline double get_cost(int i, int j) const {
            assert(i >= get_vertices_begin() && i < get_vertices_end());
            assert(j >= get_vertices_begin() && j < get_vertices_end());

            const double sqrt = std::sqrt((xcoords[i] - xcoords[j]) * (xcoords[i] - xcoords[j]) +
                                          (ycoords[i] - ycoords[j]) * (ycoords[i] - ycoords[j]));
            assert(static_cast<int>(std::round(sqrt)) == static_cast<int>(fastround(sqrt)));

            return fastround(sqrt);
        }

        // Returns the demand of vertex `i`. The demand is 0 for the depot.
        inline int get_demand(int i) const {
            return demands[i];
        };

        // Returns the x coordinate of vertex `i`.
        inline double get_x_coordinate(int i) const {
            return xcoords[i];
        };

        // Returns the y coordinate of vertex `i`.
        inline double get_y_coordinate(int i) const {
            return ycoords[i];
        };

        // Returns an array of vertices sorted according to increasing cost from `i`. This array always includes `i` in the first
        // position. The total number of elements in the array is defined by the `num_neighbors` parameter in the `Instance` constructor.
        inline const std::vector<int>& get_neighbors_of(int i) const {
            return neighbors[i];
        };

        void transform_to_unit_demands() {
            // Calcular demanda total original
            int total_demand = 0;
            for (int i = get_customers_begin(); i < get_customers_end(); ++i) {
                total_demand += demands[i];
            }
            
            // Calcular número esperado de rutas con demandas originales
            int expected_routes = (total_demand + vehicle_capacity - 1) / vehicle_capacity;
            
            int n_customers = get_customers_num();
            
            std::cout << "\n=== TRANSFORMING TO UNIT DEMANDS ===\n";
            std::cout << "Original configuration:\n";
            std::cout << "  - Customers: " << n_customers << "\n";
            std::cout << "  - Total demand: " << total_demand << "\n";
            std::cout << "  - Vehicle capacity: " << vehicle_capacity << "\n";
            std::cout << "  - Expected routes: " << expected_routes << "\n";
            
            // Transformar todas las demandas a 1
            for (int i = get_customers_begin(); i < get_customers_end(); ++i) {
                demands[i] = 1;
            }
            
            // Ajustar capacidad para mantener número aproximado de rutas
            // Nueva capacidad = ⌈n_customers / expected_routes⌉
            vehicle_capacity = (n_customers + expected_routes - 1) / expected_routes;
            
            std::cout << "\nNew configuration:\n";
            std::cout << "  - Customers: " << n_customers << "\n";
            std::cout << "  - Total demand: " << n_customers << " (all demands = 1)\n";
            std::cout << "  - Vehicle capacity: " << vehicle_capacity << "\n";
            std::cout << "  - Expected routes: " << (n_customers + vehicle_capacity - 1) / vehicle_capacity << "\n";
            std::cout << "====================================\n\n";
        }

    private:
        Instance(const Parser::Data& data, int neighbors_num);

        // Maximum vehicle capacity.
        int vehicle_capacity;

        // Vertices x coordinates.
        std::vector<double> xcoords;

        // Vertices y coordinates.
        std::vector<double> ycoords;

        // Vertices demands.
        std::vector<int> demands;

        // Neighbors for each vertex.
        std::vector<std::vector<int>> neighbors;
    };

}  // namespace cobra

#endif