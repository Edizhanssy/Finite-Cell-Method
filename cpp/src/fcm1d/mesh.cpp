#include "fcm/fcm1d/mesh.hpp"

namespace fcm {

    Mesh build_mesh(const Config& cfg) {
        Mesh m;
        const int ne         = cfg.n_elements;
        const int n_internal = cfg.p - 1;

        m.node_local.resize(static_cast<std::size_t>(ne) + 1);
        for (int i = 0; i <= ne; ++i)
            m.node_local[static_cast<std::size_t>(i)] = -1.0 + 2.0 * i / ne;

        const double dx = cfg.L / ne;
        m.elements.reserve(static_cast<std::size_t>(ne));
        for (int i = 0; i < ne; ++i)
            m.elements.push_back(Element1D{
                i,
                {m.node_local[static_cast<std::size_t>(i)],
                 m.node_local[static_cast<std::size_t>(i) + 1]},
                {i * dx, (i + 1) * dx}});

        // DOF numbering
        int counter = 0;
        m.node_dof.assign(static_cast<std::size_t>(ne) + 1, 0);
        for (int i = 0; i <= ne; ++i)
            m.node_dof[static_cast<std::size_t>(i)] = ++counter;

        m.edge_dof.assign(static_cast<std::size_t>(ne), std::vector<int>(n_internal, 0));
        for (int j = 0; j < n_internal; ++j)
            for (int e = 0; e < ne; ++e)
                m.edge_dof[static_cast<std::size_t>(e)][static_cast<std::size_t>(j)] = ++counter;
        m.n_dof = counter;

        m.ltog.assign(static_cast<std::size_t>(ne), std::vector<int>(cfg.n_modes(), 0));
        for (int e = 0; e < ne; ++e) {
            auto& L = m.ltog[static_cast<std::size_t>(e)];
            L[0] = m.node_dof[static_cast<std::size_t>(e)] - 1;
            L[1] = m.node_dof[static_cast<std::size_t>(e) + 1] - 1;
            for (int j = 0; j < n_internal; ++j)
                L[static_cast<std::size_t>(j) + 2] =
                    m.edge_dof[static_cast<std::size_t>(e)][static_cast<std::size_t>(j)] - 1;
        }
        return m;
    }

}