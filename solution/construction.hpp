#ifndef _FILO2_CONSTRUCTION_HPP_
#define _FILO2_CONSTRUCTION_HPP_

#include <vector>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <cstring>
#include <unistd.h>
#include <fcntl.h>

#include "../instance/Instance.hpp"
#include "../solution/Solution.hpp"

extern "C" {
    // Suprimir warnings de redefinición de palabras clave
    #if defined(__clang__)
        #pragma clang diagnostic push
        #pragma clang diagnostic ignored "-Wkeyword-macro"
    #elif defined(__GNUC__)
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wpragmas"
    #endif
    
    #define new _new_var_
    #define class _class_var_
    
    #ifdef NDEBUG
    #undef NDEBUG
    #endif
    
    #include "config.h"
    #include "concorde.h"
    
    #undef new
    #undef class
    
    #if defined(__clang__)
        #pragma clang diagnostic pop
    #elif defined(__GNUC__)
        #pragma GCC diagnostic pop
    #endif
    
    #include "machdefs.h"
    #include "util.h"
    #include "kdtree.h"
    #include "macrorus.h"
    #include "linkern.h"
}

#define ERROR (1e-6)

namespace cobra {
    
    // Clase para silenciar completamente stdout/stderr de Concorde
    class ConcordeSilencer {
    private:
        int saved_stdout;
        int saved_stderr;
        int dev_null;
        bool active;
        
    public:
        ConcordeSilencer() : active(true) {
            if (active) {
                // Guardar descriptores originales
                saved_stdout = dup(STDOUT_FILENO);
                saved_stderr = dup(STDERR_FILENO);
                
                // Abrir /dev/null
                dev_null = open("/dev/null", O_WRONLY);
                
                // Redirigir stdout y stderr a /dev/null
                dup2(dev_null, STDOUT_FILENO);
                dup2(dev_null, STDERR_FILENO);
            }
        }
        
        ~ConcordeSilencer() {
            if (active) {
                // Restaurar descriptores originales
                dup2(saved_stdout, STDOUT_FILENO);
                dup2(saved_stderr, STDERR_FILENO);
                
                // Cerrar descriptores temporales
                close(saved_stdout);
                close(saved_stderr);
                close(dev_null);
            }
        }
    };
    
    // Estructura para coordenadas polares (similar a la del tutor)
    struct PolarCoord {
        int customer_id;
        double rho;        
        double phi;        
        int demand;
    };

    // Variables globales para Concorde (como en el código del tutor)
    static int norm = CC_EUCLIDEAN;
    static int seed = 0;
    static int quadtry = 2;
    static int run_silently = 1;
    static int kick_type = CC_LK_WALK_KICK;
    static int current_in_repeater = 5; // MEDIUM_IN_REPEATER
    static double time_bound = -1.0;
    static double length_bound = -1.0;
    static char *saveit_name = (char *) NULL;
    
    // Implementación fiel del callLinkern del tutor
    inline double callLinkern(int n, int init, const Instance &instance, std::vector<int> &assignment) {
        int depot_id = instance.get_depot();
        
        // Casos especiales exactamente como el tutor
        if (n == 2) {
            double rho = instance.get_cost(depot_id, assignment[init]);
            return 2.0 * rho;
        }
        
        if (n == 3) {
            double rho1 = instance.get_cost(depot_id, assignment[init]);
            double rho2 = instance.get_cost(depot_id, assignment[init+1]);
            double dist = instance.get_cost(assignment[init], assignment[init+1]);
            return rho1 + rho2 + dist;
        }
        
        // Silenciar la salida de Concorde
        ConcordeSilencer silencer;
        
        // Para n >= 4, usar Concorde exactamente como el tutor
        // Declarar TODAS las variables al principio para evitar problemas con goto
        int ncount = n;
        double val = 0.0;
        int tempcount = 0;
        int *templist = nullptr;
        int *incycle = nullptr;
        int *outcycle = nullptr;
        CCdatagroup dat;
        CCrandstate rstate;
        bool kdtree_built = false;
        CCkdtree localkt;
        int in_repeater;
        int i;
        
        CCutil_init_datagroup(&dat);
        
        seed = (int) CCutil_real_zeit();
        norm = CC_EUCLIDEAN;
        CCutil_dat_setnorm(&dat, norm);
        CCutil_sprand(seed, &rstate);
        
        dat.x = CC_SAFE_MALLOC(ncount, double);
        dat.y = CC_SAFE_MALLOC(ncount, double);
        
        if (!dat.x || !dat.y) {
            goto CLEANUP;
        }
        
        // Copiar coordenadas exactamente como el tutor
        for (i = 0; i < n; ++i) {
            int cust_id = assignment[init + i];
            dat.x[i] = instance.get_x_coordinate(cust_id);
            dat.y[i] = instance.get_y_coordinate(cust_id);
        }
        
        CCutil_dat_getnorm(&dat, &norm);
        
        in_repeater = std::min(ncount, current_in_repeater);
        
        incycle = CC_SAFE_MALLOC(ncount, int);
        if (!incycle) {
            goto CLEANUP;
        }
        
        // Intentar construir kdtree y qboruvka tour
        if (CCkdtree_build(&localkt, ncount, &dat, (double *) NULL, &rstate) == 0) {
            kdtree_built = true;
            
            if (CCkdtree_quadrant_k_nearest(&localkt, ncount, quadtry, &dat, (double *) NULL, 1, 
                                             &tempcount, &templist, run_silently, &rstate) == 0) {
                
                // Intentar qboruvka, si falla usar tour trivial
                if (CCkdtree_qboruvka_tour(&localkt, ncount, &dat, incycle, &val, &rstate) != 0) {
                    // qboruvka falló, usar tour trivial: 0, 1, 2, ..., n-1
                    for (i = 0; i < ncount; ++i) {
                        incycle[i] = i;
                    }
                }
            } else {
                // quadrant_k_nearest falló, usar tour trivial
                for (i = 0; i < ncount; ++i) {
                    incycle[i] = i;
                }
            }
        } else {
            // kdtree_build falló, usar tour trivial
            for (i = 0; i < ncount; ++i) {
                incycle[i] = i;
            }
        }
        
        // Liberar kdtree solo si se construyó
        if (kdtree_built) {
            CCkdtree_free(&localkt);
        }
        
        outcycle = CC_SAFE_MALLOC(ncount, int);
        if (!outcycle) {
            goto CLEANUP;
        }
        
        // CClinkern_tour exactamente con los mismos parámetros que el tutor
        if (CClinkern_tour(ncount, &dat, tempcount, templist, 1000000000, in_repeater, 
                          incycle, outcycle, &val, run_silently, time_bound, length_bound, 
                          saveit_name, kick_type, &rstate) != 0) {
            // LK falló, usar el tour de entrada como salida
            for (i = 0; i < ncount; ++i) {
                outcycle[i] = incycle[i];
            }
            // Calcular distancia manualmente
            val = 0.0;
            for (i = 0; i < ncount; ++i) {
                int from = assignment[init + outcycle[i]];
                int to = assignment[init + outcycle[(i + 1) % ncount]];
                val += instance.get_cost(from, to);
            }
        }
        
        // Transformación de índices exactamente como el tutor (líneas 345-351)
        // Primero convertir outcycle de índices locales a IDs de clientes
        for (i = 0; i < n; ++i) {
            outcycle[i] = assignment[init + outcycle[i]];
        }
        
        // Luego copiar de vuelta al assignment
        for (i = 0; i < n; ++i) {
            assignment[init + i] = outcycle[i];
        }
        
    CLEANUP:
        CC_IFFREE(templist, int);
        CC_IFFREE(incycle, int);
        CC_IFFREE(outcycle, int);
        CCutil_freedatagroup(&dat);
        
        return val;
    }
    
    // Implementación fiel del packInTrucks del tutor
    inline void packInTrucks(const Instance &instance, 
                             const std::vector<PolarCoord> &toPack, 
                             std::vector<int> &assignment) {
        int remainingCapacity = instance.get_vehicle_capacity();
        int depot_id = instance.get_depot();
        
        for (const auto &pc : toPack) {
            if (remainingCapacity < pc.demand) {
                // Agregar depot + separador 0 exactamente como el tutor
                assignment.push_back(depot_id);
                assignment.push_back(0);
                remainingCapacity = instance.get_vehicle_capacity();
            }
            assert(remainingCapacity >= pc.demand);
            assignment.push_back(pc.customer_id);
            remainingCapacity -= pc.demand;
            assert(remainingCapacity >= 0);
        }
        
        // Agregar depot + separador 0 al final exactamente como el tutor
        assignment.push_back(depot_id);
        assignment.push_back(0);
    }
    
    // Implementación fiel del computeRoutes del tutor
    inline double computeRoutes(const Instance &instance, std::vector<int> &assignment) {
        int pointerToInit = 0;
        int n = 0;
        double totalDistance = 0;
        int demand = 0;
        
        #ifdef VERBOSE
        std::cout << "Computing routes from assignment..." << std::endl;
        #endif
        
        for (unsigned int i = 0; i < assignment.size(); ++i) {
            // El separador 0 marca el fin de una ruta (exactamente como el tutor)
            if (assignment[i] == 0) {
                if (n == 0) {
                    // Ruta vacía, saltarla
                    pointerToInit = i + 1;
                    continue;
                } else if (n == 1) {
                    // Ruta con un solo cliente: depot -> customer -> depot
                    int depot_id = instance.get_depot();
                    int customer = assignment[pointerToInit];
                    totalDistance += 2.0 * instance.get_cost(depot_id, customer);
                    pointerToInit = i + 1;
                    n = 0;
                    demand = 0;
                } else {
                    // Ruta con 2 o más clientes
                    if (demand > instance.get_vehicle_capacity()) {
                        fprintf(stderr, "assignment with more demand than capacity\n");
                        exit(10);
                    }
                    
                    #ifdef VERBOSE
                    std::cout << "Calling LK for route with " << n << " customers" << std::endl;
                    #endif
                    
                    totalDistance += callLinkern(n, pointerToInit, instance, assignment);
                    pointerToInit = i + 1;
                    n = 0;
                    demand = 0;
                }
            } else {
                ++n;
                demand += instance.get_demand(assignment[i]);
            }
        }
        
        if (totalDistance == 0 + ERROR) {
            exit(10);
        }
        
        #ifdef VERBOSE
        std::cout << "Total distance computed: " << totalDistance << std::endl;
        #endif
        
        return totalDistance;
    }
    
    // Función para construir la solución en FILO2 a partir del assignment optimizado
    inline void buildSolutionFromAssignment(const Instance &instance, 
                                            const std::vector<int> &assignment, 
                                            Solution &solution) {
        solution.reset();
        int depot_id = instance.get_depot();
        std::vector<int> currentRoute;
        
        #ifdef VERBOSE
        std::cout << "Building solution from assignment of size " << assignment.size() << std::endl;
        #endif
        
        for (unsigned int i = 0; i < assignment.size(); ++i) {
            if (assignment[i] == 0) {
                // Fin de ruta: construir la ruta en FILO2
                if (!currentRoute.empty()) {
                    #ifdef VERBOSE
                    std::cout << "Building route with " << currentRoute.size() << " customers" << std::endl;
                    #endif
                    
                    // Verificar que todos los clientes sean válidos
                    for (size_t k = 0; k < currentRoute.size(); ++k) {
                        if (currentRoute[k] < instance.get_customers_begin() || 
                            currentRoute[k] >= instance.get_customers_end()) {
                            std::cerr << "ERROR: Invalid customer ID " << currentRoute[k] << std::endl;
                            exit(1);
                        }
                    }
                    
                    // Crear ruta con el primer cliente
                    solution.build_one_customer_route</*record_action=*/false>(currentRoute[0]);
                    int route_idx = solution.get_route_index(currentRoute[0]);
                    
                    // Insertar el resto de clientes en orden
                    for (size_t k = 1; k < currentRoute.size(); ++k) {
                        solution.insert_vertex_before</*record_action=*/false>(
                            route_idx, depot_id, currentRoute[k]);
                    }
                    
                    currentRoute.clear();
                }
            } else if (assignment[i] != depot_id) {
                // Es un cliente, agregarlo a la ruta actual
                currentRoute.push_back(assignment[i]);
            }
            // Si es depot_id pero no es 0, lo ignoramos (es el depot entre rutas)
        }
        
        #ifdef VERBOSE
        std::cout << "Solution built with " << solution.get_routes_num() << " routes" << std::endl;
        #endif
        
        assert(solution.is_feasible());
    }
    
    // Implementación fiel del sweep del tutor
    inline void sweep(const Instance &instance, Solution &solution, double th, int rotation_idx = 0) {
        #ifdef VERBOSE
        std::cout << "Starting sweep with th=" << th << ", rotation_idx=" << rotation_idx << std::endl;
        std::cout << "Customers: " << instance.get_customers_num() << ", Capacity: " << instance.get_vehicle_capacity() << std::endl;
        #endif
        
        solution.reset();
        
        // Calcular coordenadas polares para todos los clientes
        std::vector<PolarCoord> polar_coords;
        polar_coords.reserve(instance.get_customers_num());
        
        const auto depot = instance.get_depot();
        const auto depot_x = instance.get_x_coordinate(depot);
        const auto depot_y = instance.get_y_coordinate(depot);
        
        #ifdef VERBOSE
        std::cout << "Depot: " << depot << " at (" << depot_x << ", " << depot_y << ")" << std::endl;
        #endif
        
        // Calcular rhoMax mientras construimos las coordenadas polares
        double rhoMax = 0.0;
        
        for (auto i = instance.get_customers_begin(); i < instance.get_customers_end(); i++) {
            const auto dx = instance.get_x_coordinate(i) - depot_x;
            const auto dy = instance.get_y_coordinate(i) - depot_y;
            
            PolarCoord pc;
            pc.customer_id = i;
            pc.rho = std::sqrt(dx * dx + dy * dy);
            pc.phi = std::atan2(dy, dx);
            pc.demand = instance.get_demand(i);
            
            polar_coords.push_back(pc);
            
            if (pc.rho > rhoMax) {
                rhoMax = pc.rho;
            }
        }
        
        #ifdef VERBOSE
        std::cout << "Computed polar coords for " << polar_coords.size() << " customers, rhoMax=" << rhoMax << std::endl;
        #endif
        
        // Ordenar por ángulo phi (como orderedByPhi del tutor)
        std::sort(polar_coords.begin(), polar_coords.end(), 
                  [](const PolarCoord& a, const PolarCoord& b) { 
                      return a.phi < b.phi; 
                  });
        
        // NOTA: rotation_idx está presente para compatibilidad con main.cpp
        // pero NO se usa en esta implementación fiel al tutor
        // La rotación se implementará en una versión futura
        (void)rotation_idx; // Suprimir warning de variable no usada
        
        // Calcular threshold EXACTAMENTE como el tutor (línea 200)
        double filterValue = th * rhoMax + ERROR;
        
        #ifdef VERBOSE
        std::cout << "Filter value: " << filterValue << std::endl;
        #endif
        
        // Separar en inner y outer según el threshold
        std::vector<PolarCoord> inner, outer;
        for (const auto& pc : polar_coords) {
            if (pc.rho < filterValue) {
                inner.push_back(pc);
            } else {
                outer.push_back(pc);
            }
        }
        
        #ifdef VERBOSE
        std::cout << "Inner: " << inner.size() << ", Outer: " << outer.size() << std::endl;
        #endif
        
        // Crear assignment exactamente como el tutor
        std::vector<int> assignment;
        if (!inner.empty()) {
            packInTrucks(instance, inner, assignment);
        }
        if (!outer.empty()) {
            packInTrucks(instance, outer, assignment);
        }
        
        #ifdef VERBOSE
        std::cout << "Assignment size: " << assignment.size() << std::endl;
        std::cout << "Calling computeRoutes..." << std::endl;
        #endif
        
        // Optimizar rutas con LK
        computeRoutes(instance, assignment);
        
        #ifdef VERBOSE
        std::cout << "Routes computed, building solution..." << std::endl;
        #endif
        
        // Construir la solución en FILO2
        buildSolutionFromAssignment(instance, assignment, solution);
        
        #ifdef VERBOSE
        std::cout << "Solution built successfully" << std::endl;
        #endif
        
        assert(solution.is_feasible());
    }
    
} // namespace cobra

#endif