import gmsh
import pyvista as pv
import os

os.chdir(os.path.dirname(__file__))

# Initialize the gmsh session
gmsh.initialize()

# Load your .msh file which contains the original model and mesh data
gmsh.open('transformed_atm.msh')

# Tell Gmsh to re-generate the 2D mesh
gmsh.model.mesh.generate(2)

# Save the mesh to a temporary .msh file
temp_mesh_path = 'temp_mesh.msh'
gmsh.write(temp_mesh_path)

# Finalize the gmsh session
gmsh.finalize()

# Now, load the temporary mesh file with PyVista
mesh = pv.read(temp_mesh_path)

# Convert mesh to VTU and save it
vtu_mesh_path = 'transformed_atm.vtu'
mesh.save(vtu_mesh_path)

# Optionally, delete the temporary .msh file
os.remove(temp_mesh_path)