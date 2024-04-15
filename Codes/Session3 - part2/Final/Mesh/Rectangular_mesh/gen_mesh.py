import gmsh
import os

# Set the current directory to the directory containing the script
script_directory = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_directory)

def create_structured_rectangle(width, height):
    gmsh.initialize()
    gmsh.model.add("structured_rectangle")

    # Define the rectangle
    p1 = gmsh.model.geo.addPoint(0, 0, 0)
    p2 = gmsh.model.geo.addPoint(width, 0, 0)
    p3 = gmsh.model.geo.addPoint(width, height, 0)
    p4 = gmsh.model.geo.addPoint(0, height, 0)

    # Define lines
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)

    # Curve loop and plane surface
    cl = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    s = gmsh.model.geo.addPlaneSurface([cl])

    # Apply transfinite meshing with a progression factor
    numPointsLine1and3 = 11 # Number of divisions on lines parallel to the width
    numPointsLine2and4 = 11  # Number of divisions on lines parallel to the height

    # Set transfinite meshing on the lines with a progression
    gmsh.model.geo.mesh.setTransfiniteCurve(l1, numPointsLine1and3)
    gmsh.model.geo.mesh.setTransfiniteCurve(l2, numPointsLine2and4)
    gmsh.model.geo.mesh.setTransfiniteCurve(l3, numPointsLine1and3)
    gmsh.model.geo.mesh.setTransfiniteCurve(l4, numPointsLine2and4)

    # Set the entire surface as transfinite, forcing a structured grid
    gmsh.model.geo.mesh.setTransfiniteSurface(s, "Right", [l1, l2, l3, l4])

    # Optionally, recombine the mesh to form quadrilaterals
    gmsh.model.geo.mesh.setRecombine(2, s)

    gmsh.model.geo.synchronize()

    # Generate the mesh
    gmsh.model.mesh.generate(2)
   
    # Specify the MSH version to export in MSH 2.2 format # CRUCIAL FOR MFEM.
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    
    # Save and view
    gmsh.write("structured_rectangle.msh")
    gmsh.fltk.run()
    
    gmsh.finalize()

# Specify the dimensions of the rectangle
width = 1
height = 1
create_structured_rectangle(width, height)



