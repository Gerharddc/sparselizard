#include "sparselizard.h"
#include "gmsh.h"

using namespace sl;

int main(void)
{
    // The domain regions as defined in the mesh file
    int substrate = 57, all_copper = 58, in_face = 59, out_face = 60;

    std::cout << "Start mesh load\n";
    gmsh::initialize();
    gmsh::open("SimplePCB.msh");
    // Load mesh available in GMSH API.
    // Set verbosity to 2 to print the physical regions info.
    mesh mymesh("gmsh:api", 2);
    gmsh::finalize();
    std::cout << "Finish mesh load\n";

    // Field V is the electric potential
    field V("h1");

    V.setorder(all_copper, 2);

    // Define boundary conditions
    V.setconstraint(in_face, 5.0);  // Input voltage is 5V
    V.setconstraint(out_face, 0.0); // Output voltage is 5V

    // Define Material Properties
    const double sigma_copper = 5.9e7; // Conductivity in S/m

    // Set Up the Formulation
    formulation voltage_drop;
    voltage_drop += integral(all_copper, -sigma_copper * grad(dof(V)) * grad(tf(V)));

    // Generate, solve and transfer the solution to field u:
    std::cout << "Starting solve\n";
    voltage_drop.solve();
    std::cout << "Finished solve\n";

    // Post-Processing
    V.write(all_copper, "voltage_drop.vtk", 2);
}
