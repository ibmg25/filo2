#ifndef _FILO2_CONSTRUCTION_HPP_
#define _FILO2_CONSTRUCTION_HPP_

#include "../instance/Instance.hpp"
#include "Solution.hpp"

namespace cobra {

    // Implementación del Sweep con descomposición radial y optimización Linkern (Concorde)
    // Se prueba un rango de thresholds para separar nodos internos/externos.
    void sweep_concorde_heuristic(const Instance &instance, Solution &solution, int seed);

}

#endif