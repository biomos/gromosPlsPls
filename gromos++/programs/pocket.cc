/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file pocket.cc
 * Calculates the pocket around the input specified coordinate point. 
 */

/**
 * @page programs Program Documentation
 *
 * @anchor pocket
 * @section Calculates the pocket around the input specified coordinate point
 * @author @ref tk
 * @date 10-04-25
 *
 * This program calculates the binding site vectors spreading circularly from the 
 * chosen ancor/central point. 
 * 
 * Prerequirements:
 * The snapshots of the trajectory you are investigating have to be aligned.
 * Once you know where your pocket is, find the coordinate that is preferably the 
 * center of your cavity. The center point is read in using vector specifier, so
 * it can be read in in format of cartesian coordinates, polar coordinates, or as 
 * a selected atom of the system.
 *
 * How it works:
 * Before the vectors are calculated, a sphere with an input-defined radius (nm) and number 
 * of vectors will be generated. This sphere is then translated to have its center in 
 * the center of your cavity that you defined. In case you are interested in half of
 * sphere only, use the flag @hemisphere to keep only the sphere part with z>=0.
 * For every snapshot of the trajectory, new vector lengths will be generated so that 
 * every vector stops once it hits an atom of the protein. Lastly, vectors will also 
 * be asigned the partial charge of the atom they hit. If no atom was encountered,
 * the vector remains at its full length and the charge is set to 0.
 *
 * Default output files include:
 * - lengths.txt: inludes the final lengths for all the vectors, one snapshot per line
                  with first element of each line containing the snapshot order number (SNAPSHOT_0, SNAPSHOT_1...)
 * - charges.txt: inludes the final charges for all the vectors, one snapshot per line
                  with first element of each line containing the snapshot order number (SNAPSHOT_0, SNAPSHOT_1...)

 * Optional output files include:
 * - volume.txt: shows the volume of the pocket
 * - area.txt: shows the area of the pocket
 * - the volume of the pocket for every snapshot
 * - the area of the pocket for every snapshot
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gathermethod&gt;] </td></tr>
 * <tr><td> \@center</td><td>&lt;pocket center&gt; </td></tr>
 * <tr><td> \@protein</td><td>&lt;protein atoms to consider for vector generation&gt; </td></tr>
 * <tr><td> [\@reject</td><td>&lt;protein atoms to discard for vector generation] </td></tr>
 * <tr><td> \@radius</td><td>&lt;max. length of the generated binding site vectors (nm) &gt; </td></tr>
 * <tr><td> \@vec_number_factor</td><td>&lt;set factor that decides the number of binding site vectors&gt; </td></tr>
 
 * <tr><td> [\@radH</td><td>&lt;radius to be used for hydrogen atoms; 0.11 nm by default] </td></tr>
 * <tr><td> [\@hemisphere</td><td>&lt;keep only the z>0 hemisphere initial vectors] </td></tr>
 * <tr><td> [\@volume_and_area</td><td>&lt;compute enclosed volume and surface area] </td></tr>
 * <tr><td> [\@final_vector_coords</td><td>&lt;coordinates of the truncated vectors] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  pocket
    @topo               topology.top
    @pbc                r cog
    @center             'cart(0,0,0)' OR 'atom(1:4875)'
    @protein            1:a
    @reject             1:4875-4921 
    @radius             2
    @vec_number_factor  7
    @radH               0.0
    @hemisphere
    @volume_and_area
    @final_vector_coords
    @traj               traj.trc
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <vector>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <cmath>
#include <set>
#include <algorithm>
#include <fstream>
#include <limits>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/LJException.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/utils/Value.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/groTime.h"
#include "../src/utils/AtomicRadii.h"
#include "../src/utils/VectorSpecifier.h"

using namespace fit;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace std;
using namespace utils;

using Vec3 = std::array<double, 3>;
constexpr double phi = (1.0 + std::sqrt(5.0)) / 2.0;

Vec3 normalize(const Vec3& v) {
    double norm = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    return {v[0]/norm, v[1]/norm, v[2]/norm};
}

Vec3 linear_combination(double u, double v, double w, const Vec3& v1, const Vec3& v2, const Vec3& v3) {
    return {
        u*v1[0] + v*v2[0] + w*v3[0],
        u*v1[1] + v*v2[1] + w*v3[1],
        u*v1[2] + v*v2[2] + w*v3[2]
    };
}

Vec3 round_point(const Vec3& p, int decimals = 8) {
    double factor = std::pow(10, decimals);
    return {
        std::round(p[0]*factor)/factor,
        std::round(p[1]*factor)/factor,
        std::round(p[2]*factor)/factor
    };
}

struct Vec3Compare {
    bool operator()(const Vec3& a, const Vec3& b) const {
        if (a[0] != b[0]) return a[0] < b[0];
        if (a[1] != b[1]) return a[1] < b[1];
        return a[2] < b[2];
    }
};

std::vector<Vec3> icosahedron_vertices() {
    std::vector<Vec3> vertices = {
        {-1, phi, 0}, {1, phi, 0}, {-1, -phi, 0}, {1, -phi, 0},
        {0, -1, phi}, {0, 1, phi}, {0, -1, -phi}, {0, 1, -phi},
        {phi, 0, -1}, {phi, 0, 1}, {-phi, 0, -1}, {-phi, 0, 1}
    };
    for (auto& v : vertices)
        v = normalize(v);
    return vertices;
}

std::vector<std::array<int, 3>> icosahedron_faces() {
    return {
        {0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11},
        {1, 5, 9}, {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
        {3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8}, {3, 8, 9},
        {4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1}
    };
}

std::vector<Vec3> subdivide_triangle(const Vec3& v1, const Vec3& v2, const Vec3& v3, int n) {
    std::vector<Vec3> points;
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= n - i; ++j) {
            double u = double(i) / n;
            double v = double(j) / n;
            double w = 1.0 - u - v;
            Vec3 p = linear_combination(u, v, w, v1, v2, v3); 
            p = normalize(p);
            points.push_back(p);
        }
    }
    return points;
}

// add includes at top if not already:
// #include <map>
// #include <tuple>

struct SphereMesh {
    std::vector<Vec3> vertices;
    std::vector<std::array<int, 3>> triangles;
};

// helper to compute triangle area & tetra volume (origin)
double triangle_area(const Vec3& v1, const Vec3& v2, const Vec3& v3) {
    Vec3 a = {v2[0]-v1[0], v2[1]-v1[1], v2[2]-v1[2]};
    Vec3 b = {v3[0]-v1[0], v3[1]-v1[1], v3[2]-v1[2]};
    Vec3 cr = {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
    double norm = std::sqrt(cr[0]*cr[0] + cr[1]*cr[1] + cr[2]*cr[2]);
    return 0.5 * norm;
}

double pyramid_volume(const Vec3& v1, const Vec3& v2, const Vec3& v3) {
    // Volume of tetrahedron from origin to triangle (v1,v2,v3): |dot(v1, cross(v2, v3))|/6
    double vol = (v1[0]*(v2[1]*v3[2] - v2[2]*v3[1])
                - v1[1]*(v2[0]*v3[2] - v2[2]*v3[0])
                + v1[2]*(v2[0]*v3[1] - v2[1]*v3[0])) / 6.0;
    return std::fabs(vol);
}

// robust sphere mesh generator: returns vertices (radius-scaled) and triangles (indices)
SphereMesh generate_sphere_mesh(double radius, int subdivisions, const Vec3& /*center*/) {
    // subdivisions = n (number subdivisions per edge)
    if (subdivisions < 1) subdivisions = 1;

    auto ico_vertices = icosahedron_vertices();
    auto faces = icosahedron_faces();

    // map rounded unit direction -> global vertex index
    std::map<Vec3, int, Vec3Compare> vertex_map;
    std::vector<Vec3> unique_vertices;
    std::vector<std::array<int, 3>> triangles;
    // for each icosa face create a triangular grid of (n+1)(n+2)/2 points
    int n = subdivisions;
    for (const auto &face : faces) {
        Vec3 v1 = ico_vertices[face[0]];
        Vec3 v2 = ico_vertices[face[1]];
        Vec3 v3 = ico_vertices[face[2]];

        // grid to store indices for this face: idxGrid[i][j] where 0 <= i <= n, 0 <= j <= n-i
        std::vector<std::vector<int>> idxGrid(n+1);
        for (int i = 0; i <= n; ++i) {
            idxGrid[i].resize(n + 1 - i, -1);
        }

        // create points & map to unique index
        for (int i = 0; i <= n; ++i) {
            for (int j = 0; j <= n - i; ++j) {
                double u = double(i) / double(n);
                double v = double(j) / double(n);
                double w = 1.0 - u - v;
                Vec3 p = linear_combination(u, v, w, v1, v2, v3);
                p = normalize(p);
                Vec3 rp = round_point(p, 8); // canonical key

                auto it = vertex_map.find(rp);
                int idx;
                if (it == vertex_map.end()) {
                    idx = unique_vertices.size();
                    vertex_map[rp] = idx;
                    // scale to radius and keep centered at origin (we translate system elsewhere)
                    unique_vertices.push_back({rp[0] * radius, rp[1] * radius, rp[2] * radius});
                } else {
                    idx = it->second;
                }
                idxGrid[i][j] = idx;
            }
        }

        // build triangles within the grid
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n - i; ++j) {
                int a = idxGrid[i][j];
                int b = idxGrid[i+1][j];
                int c = idxGrid[i][j+1];
                // triangle (a,b,c)
                triangles.push_back({a,b,c});

                // second triangle in the square: (b, idxGrid[i+1][j+1], c) if exists
                if (j + i + 1 <= n - 1) {
                    int d = idxGrid[i+1][j+1];
                    triangles.push_back({b,d,c});
                }
            }
        }
    }

    SphereMesh mesh;
    mesh.vertices = std::move(unique_vertices);
    mesh.triangles = std::move(triangles);
    return mesh;
}
// validation helper: compute mesh area/volume and print relative error to analytic sphere
void validate_mesh(const SphereMesh &mesh, double radius) {
    double area_sum = 0.0;
    double vol_sum = 0.0;

    // sanity checks
    size_t nv = mesh.vertices.size();
    size_t nt = mesh.triangles.size();
    std::cerr << "Validation: mesh has " << nv << " vertices and " << nt << " triangles\n";

    // check triangle indices bounds
    for (const auto &t : mesh.triangles) {
        if ((size_t)t[0] >= nv || (size_t)t[1] >= nv || (size_t)t[2] >= nv) {
            std::cerr << "ERROR: triangle index out of range\n";
            return;
        }
    }

    for (const auto &tri : mesh.triangles) {
        const Vec3 &v1 = mesh.vertices[tri[0]];
        const Vec3 &v2 = mesh.vertices[tri[1]];
        const Vec3 &v3 = mesh.vertices[tri[2]];
        area_sum += triangle_area(v1, v2, v3);
        vol_sum += pyramid_volume(v1, v2, v3);
    }

    
    double analytic_area = 4.0 * M_PI * radius * radius;
    double analytic_vol = 4.0/3.0 * M_PI * radius * radius * radius;
    std::cerr << std::fixed << std::setprecision(6);
    std::cerr << "Mesh area of generated initial hemisphere= " << area_sum << "\n" << "Mesh area of generated initial full sphere = " << analytic_area << "\n";
    std::cerr << "Mesh volume of generated initial hemisphere  = " << vol_sum << "\n" << "Mesh volume of generated initial full sphere = " << analytic_vol << "\n";
}

int main(int argc, char **argv) {
    Argument_List knowns;
    knowns << "topo" << "pbc" << "center" << "protein"
           << "reject" << "radius" << "vec_number_factor"
           << "radH" << "hemisphere" << "volume_and_area"  
           << "final_vector_coords" << "traj";

    string usage = argv[0];
    usage += "\n\t@topo                 <molecular topology file>\n";
    usage += "\t@pbc                    <boundary type> [<gathermethod>]\n";
    usage += "\t@center                 <pocket center>\n";
    usage += "\t@protein                <protein atoms to consider for vector generation>\n";
    usage += "\t@[reject                <protein atoms to discard for vector generation>]\n";
    usage += "\t@radius                 <max. length of the generated binding site vectors (nm)>\n";
    usage += "\t@vec_number_factor      <factor that decides the number of binding site vectors>\n";
    usage += "\t[@radH                  <radius to be used for hydrogen atoms>]\n";
    usage += "\t[@hemisphere            <keep only the z>0 hemisphere initial vectors>]\n";
    usage += "\t[@volume_and_area       <compute enclosed volume and surface area>]\n";
    usage += "\t[final_vector_coords    <coordinates of the truncated vectors>]\n";
    usage += "\t@traj                   <trajectory files>\n";

    try {
        Arguments args(argc, argv, knowns, usage);

        InTopology it(args["topo"]);
        System sys(it.system());

        Boundary *pbc = BoundaryParser::boundary(sys, args);
        Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, sys, args);

        bool use_hemisphere = false;
        if (args.count("hemisphere") >= 0)
            use_hemisphere = true;

        bool volume_and_area = false;
        if (args.count("volume_and_area") >= 0)
            volume_and_area = true;
        
        bool final_vector_coords = false;
        if (args.count("final_vector_coords") >= 0)
                final_vector_coords = true;   


        utils::VectorSpecifier vs(sys, pbc);
        vs.setSpecifier(args["center"]);

        utils::AtomSpecifier protein(sys);
        {
            Arguments::const_iterator iter = args.lower_bound("protein");
            Arguments::const_iterator to = args.upper_bound("protein");
            for (; iter != to; iter++) {
                protein.addSpecifier(iter->second.c_str());
            }
        }

        // Handle reject atoms
        utils::AtomSpecifier rej(sys);
        if (args.count("reject") > 0) {
            Arguments::const_iterator iter = args.lower_bound("reject");
            Arguments::const_iterator to = args.upper_bound("reject");
            for (; iter != to; ++iter) {
                rej.addSpecifier(iter->second.c_str());
            }
        }

        // Remove rejected atoms from protein
        for (int i = 0; i < rej.size(); ++i) {
            int idx = protein.findAtom(rej.mol(i), rej.atom(i));
            if (idx != -1)
                protein.removeAtom(idx);
        }

        

        double radius = args.getValue<double>("radius", true);
        int vec_number_factor = args.getValue<int>("vec_number_factor", false, 0);
        double radH = args.getValue<double>("radH", false, 0.11);

        utils::compute_atomic_radii_vdw(sys, it.forceField());

        // Set radii for H atoms
        for (int m = 0; m < sys.numMolecules(); ++m) {
            auto &mol = sys.mol(m);
            for (int a = 0; a < mol.topology().numAtoms(); ++a) {
                auto &atomTop = mol.topology().atom(a);
                if (atomTop.radius() == 0.0)
                    atomTop.setradius(radH);
            }
        }


        gmath::Vec center = vs();

        // Generate sphere points
        SphereMesh sphere_mesh = generate_sphere_mesh(radius, vec_number_factor, {center[0], center[1], center[2]});
        std::vector<Vec3>& sphere_points = sphere_mesh.vertices;
        auto& triangles = sphere_mesh.triangles;
        std::cout << "Generated " << sphere_points.size() << " points and " << triangles.size() << " triangles." << std::endl;

        if (use_hemisphere) {
            std::vector<int> old_to_new(sphere_points.size(), -1);
            std::vector<Vec3> filtered_points;
            filtered_points.reserve(sphere_points.size());

            // Build mapping: keep only points with z ≥ 0
            for (size_t i = 0; i < sphere_points.size(); ++i) {
                if (sphere_points[i][2] >= 0.0) {
                    old_to_new[i] = static_cast<int>(filtered_points.size());
                    filtered_points.push_back(sphere_points[i]);
                }
            }

            std::vector<std::array<int, 3>> filtered_triangles;
            filtered_triangles.reserve(triangles.size());

            // Keep only triangles whose 3 vertices are in upper hemisphere
            for (const auto &tri : triangles) {
                int i1 = old_to_new[tri[0]];
                int i2 = old_to_new[tri[1]];
                int i3 = old_to_new[tri[2]];

                if (i1 != -1 && i2 != -1 && i3 != -1) {
                    filtered_triangles.push_back({i1, i2, i3});
                }
            }

            // Replace with filtered versions
            sphere_points = std::move(filtered_points);
            triangles = std::move(filtered_triangles);

            std::cout << "Hemisphere flag active: kept "
                    << sphere_points.size() << " points and "
                    << triangles.size() << " triangles above xy-plane." << std::endl;
        }




        // validate immediately (prints to stderr)
        validate_mesh(sphere_mesh, radius);



        // Output files
        std::ofstream charge_out("charges.txt");
        std::ofstream length_out("lengths.txt");
        std::ofstream volume_out;
        std::ofstream area_out;
        std::ofstream coords_out;

        if (volume_and_area) {
            volume_out.open("volume.txt");
            area_out.open("area.txt");
        }

        if (final_vector_coords)
            coords_out.open("vector_coords.txt");



        int snapshot_counter = 0;

        InG96 ic;
        for (Arguments::const_iterator iter = args.lower_bound("traj");
             iter != args.upper_bound("traj"); ++iter) {
            ic.open(iter->second);

            while (!ic.eof()) {
                ic.select("ALL");
                ic >> sys;

                gmath::Vec translation = -1.0 * center;
                for (int m = 0; m < sys.numMolecules(); ++m) {
                    Molecule &mol = sys.mol(m);
                    for (int a = 0; a < mol.numAtoms(); ++a) {
                        mol.pos(a) += translation;
                    }
                }

                std::vector<gmath::Vec> protein_coords;
                std::vector<double> radii;
                std::vector<double> charges;
                for (int i = 0; i < protein.size(); ++i) {
                    protein_coords.push_back(protein.pos(i));
                    radii.push_back(protein.radius(i));
                    charges.push_back(protein.charge(i));
                }

                std::vector<double> snapshot_distances;
                std::vector<double> snapshot_charges;
                std::vector<Vec3> final_coords;

                for (const Vec3 &vec : sphere_points) {
                    gmath::Vec V(vec[0], vec[1], vec[2]); // initial vector
                    double min_dist = std::numeric_limits<double>::infinity();
                    double t1 = std::numeric_limits<double>::infinity();
                    Vec3 final_vector = {V[0], V[1], V[2]};
                    double selected_charge = 0.0;

                    for (size_t i = 0; i < protein_coords.size(); ++i) {
                        gmath::Vec S = protein_coords[i];
                        double r = radii[i];
                        double d = S.abs();

                        if (d < 1e-6) continue; // avoid division by zero

                        double cos_alpha = S.dot(V) / (d * V.abs());
                        bool sphere_touching = false;
                        gmath::Vec temp_vec;
                        double t2 = 0.0;

                        if (S == V) {
                            t1 = d - r;
                            sphere_touching = true;
                            temp_vec = V / V.abs() * t1;  // equivalent to new_vector(V, t1)
                        }
                        else if (cos_alpha >= 0.0 && cos_alpha <= 1.0) {
                            double y = d * std::sqrt(1.0 - cos_alpha * cos_alpha);
                            if (y < r) {
                                double x = std::sqrt(r * r - y * y);
                                double t = std::abs(S.dot(V)) / V.abs();
                                t1 = t - x;
                                t2 = t + x;

                                // Limit to maximum radius
                                double used_t = (t1 < radius) ? t1 : radius;

                                temp_vec = V / V.abs() * used_t;

                                // Hemisphere enforcement (keep z ≥ 0)
                                if (temp_vec[2] >= 0.0) {
                                    sphere_touching = true;
                                } else {
                                    // For debugging you could print similar messages
                                    // std::cout << "Rejected: z<0 for vec " << i << std::endl;
                                    sphere_touching = false;
                                }
                            }
                            else if (std::abs(y - r) < 1e-6) {
                                t1 = S.dot(V) / V.abs();
                                double used_t = (t1 < radius) ? t1 : radius;
                                temp_vec = V / V.abs() * used_t;
                                sphere_touching = true;
                            }
                        }

                        if (sphere_touching && t1 < min_dist) {
                            min_dist = t1;
                            final_vector = {temp_vec[0], temp_vec[1], temp_vec[2]};
                            selected_charge = charges[i];
                        }
                    }

                    // If no intersection found, use max sphere radius
                    if (!std::isfinite(min_dist) || min_dist > radius) {
                        min_dist = radius;
                        final_vector = {V[0] / V.abs() * radius,
                                        V[1] / V.abs() * radius,
                                        V[2] / V.abs() * radius};
                    }

                    snapshot_distances.push_back(std::round(min_dist * 1000.0) / 1000.0);
                    snapshot_charges.push_back(selected_charge);
                    final_coords.push_back(final_vector);
                }

                if (volume_and_area) {
                    double total_area = 0.0;
                    double total_volume = 0.0;

                    for (const auto& tri : triangles) {
                        const Vec3& v1 = final_coords[tri[0]];
                        const Vec3& v2 = final_coords[tri[1]];
                        const Vec3& v3 = final_coords[tri[2]];


                        total_area += triangle_area(v1, v2, v3);
                        total_volume += pyramid_volume(v1, v2, v3);
                    }

                    // Write snapshot with leading zeros
                    std::string snapshot_str = std::to_string(snapshot_counter);
                    while (snapshot_str.size() < 5) snapshot_str = "0" + snapshot_str;

                    area_out << "SNAPSHOT_" << snapshot_str << " " << std::fixed << std::setprecision(6) << total_area << "\n";
                    volume_out << "SNAPSHOT_" << snapshot_str << " " << std::fixed << std::setprecision(6) << total_volume << "\n";
                }

                

                // Write to text files
                charge_out << "SNAPSHOT_" << snapshot_counter;
                length_out << "SNAPSHOT_" << snapshot_counter;
                if (final_vector_coords) 
                    coords_out << "SNAPSHOT_" << snapshot_counter << "\n";

                for (size_t i = 0; i < snapshot_charges.size(); ++i) {
                    charge_out << " " << snapshot_charges[i];
                    length_out << " " << snapshot_distances[i];

                    if (final_vector_coords) {

                    // write vector coords
                    coords_out << std::fixed << std::setprecision(6)
                               << final_coords[i][0] << ";"
                               << final_coords[i][1] << ";"
                               << final_coords[i][2] << "\n";
                    }
                }

                charge_out << "\n";
                length_out << "\n";
                
                
                snapshot_counter++;
            }
        }

        return 0;
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Unknown error occurred" << std::endl;
        return 1;
    }
}
