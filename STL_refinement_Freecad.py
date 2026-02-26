import FreeCAD, Mesh, Part
import os

# ------------------------
# USER SETTINGS
# ------------------------
input_folder  = r"D:\Compiled_NASA\VoronoiSTL"

output_folder = input_folder  # save .step in same folder (change if needed)

tolerance = 0.1   # mesh → shape tolerance

# ------------------------
# ENSURE OUTPUT FOLDER EXISTS
# ------------------------
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# ------------------------
# FIND ALL STL FILES
# ------------------------
stl_files = [f for f in os.listdir(input_folder) if f.lower().endswith(".stl")]

print(f"Found {len(stl_files)} STL files.")

# ------------------------
# LOOP THROUGH ALL STL FILES
# ------------------------
for stl_name in stl_files:
#
    stl_path = os.path.join(input_folder, stl_name)
    step_path = os.path.join(output_folder, stl_name.replace(".stl", ".step"))
#
    print("\n==========================================")
    print(f"Processing: {stl_name}")
    print("==========================================")
#
    # create new document for each conversion
    doc = FreeCAD.newDocument(stl_name.replace(".stl", ""))
#
    try:
        # ------------------------
        # LOAD STL AS MESH
        # ------------------------
        mesh = Mesh.Mesh(stl_path)
        mesh_obj = doc.addObject("Mesh::Feature", "ImportedMesh")
        mesh_obj.Mesh = mesh
        doc.recompute()
#
        # ------------------------
        # MESH → SHAPE
        # ------------------------
        shape = Part.Shape()
        shape.makeShapeFromMesh(mesh.Topology, tolerance)
#
        # Healing operations
        shape = shape.removeSplitter()
        shape.fix(tolerance, tolerance, tolerance)
#
        shape_obj = doc.addObject("Part::Feature", "MeshShape")
        shape_obj.Shape = shape
        doc.recompute()
#
        # ------------------------
        # SHAPE → SOLID
        # ------------------------
        if not shape.Shells:
            print("ERROR: No shells found → cannot convert to solid.")
            continue
#
        shell = shape.Shells[0]
        solid = Part.makeSolid(shell)
#
        # Refine
        solid = solid.removeSplitter()
#
        solid_obj = doc.addObject("Part::Feature", "Solid")
        solid_obj.Shape = solid
        doc.recompute()
#
        # ------------------------
        # VALIDATE
        # ------------------------
        valid = solid.isValid()
        print("Solid validity:", valid)
#
        if not valid:
            print("WARNING: Solid invalid — STEP export may fail.")
#
        # ------------------------
        # EXPORT STEP
        # ------------------------
        Part.export([solid_obj], step_path)
        print(f"Exported STEP: {step_path}")
#
    except Exception as e:
        print(f"FAILED on {stl_name}: {e}")
#
    finally:
        FreeCAD.closeDocument(doc.Name)

print("\n===============================")
print("BATCH CONVERSION FINISHED")
print("===============================")
