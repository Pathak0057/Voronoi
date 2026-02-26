from abaqus import *
from abaqusConstants import *
import os, re
import mesh
import math

# ============================================================
# USER SETTINGS
# ============================================================
geom_folder = r"D:\Compiled_NASA\VoronoiSTL"  # folder containing region_XXX step/stl
model_name  = "Model-1"

# Your existing CSV (must already exist)
materials_file = os.path.join(geom_folder, "materials_db.csv")

# ============================================================
# Ensure model exists
# ============================================================
if model_name not in mdb.models.keys():
    mdb.Model(name=model_name)
model = mdb.models[model_name]


# ============================================================
# Read materials CSV into dict (EXISTING FILE)
# Expected CSV format (18 columns):
# region_id,
# E1,E2,E3,nu12,nu13,nu23,G12,G13,G23,
# p10,p11,p12,
# density,conductivity,expansion,specificHeat,depvar
# ============================================================
def read_materials_csv(csv_path):
    if not os.path.isfile(csv_path):
        raise RuntimeError("materials_db.csv not found at:\n" + csv_path)

    db = {}
    with open(csv_path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue

            parts = [p.strip() for p in s.split(",")]
            if len(parts) < 18:
                raise ValueError(
                    "Bad line (expected >=18 columns).\n"
                    "Line:\n" + s + "\n\n"
                    "Expected columns:\n"
                    "region_id,E1,E2,E3,nu12,nu13,nu23,G12,G13,G23,p10,p11,p12,density,conductivity,expansion,specificHeat,depvar"
                )

            rid = parts[0]

            # Convert the rest to float; depvar handled separately
            nums = [float(x) for x in parts[1:17]]  # 16 values (9+3+4)
            depvar = int(float(parts[17]))          # depvar can be written as 25 or 25.0

            eng = tuple(nums[0:9])         # E1..G23 (9)
            umat_extra = tuple(nums[9:12]) # p10,p11,p12 (3)
            density = nums[12]
            cond    = nums[13]
            expn    = nums[14]
            cp      = nums[15]

            umat_consts = tuple(float(x) for x in (eng + umat_extra))  # FLAT tuple

            db[rid] = {
                "UMAT": umat_consts,
                "Density": float(density),
                "Conductivity": float(cond),
                "Expansion": float(expn),
                "SpecificHeat": float(cp),
                "Depvar": depvar
            }
    return db


material_database = read_materials_csv(materials_file)
print("Loaded materials from:", materials_file)
print("Number of regions in CSV:", len(material_database))


# ============================================================
# OPTIONAL: Import STEP files (if you use STEP)
# If you already imported geometry and parts exist, you can remove this block.
# ============================================================
step_files = [f for f in os.listdir(geom_folder) if f.lower().endswith((".step", ".stp"))]
print("Found STEP files:", step_files)

for step_file in step_files:
    full_path = os.path.join(geom_folder, step_file)
    part_name = os.path.splitext(step_file)[0]  # e.g., region_008

    print("\nImporting STEP:", part_name)
    mdb.openStep(full_path, scaleFromFile=OFF)

    try:
        model.PartFromGeometryFile(
            name=part_name,
            geometryFile=mdb.acis,
            dimensionality=THREE_D,
            type=DEFORMABLE_BODY,
            combine=False
        )
        print("Imported:", part_name)
    except Exception as e:
        print("WARNING: Could not import (maybe already exists):", part_name)
        print("Reason:", str(e))


# ============================================================
# CREATE UMAT MATERIALS + SECTIONS + ASSIGNMENTS
# ============================================================
print("\n=== CREATING UMAT MATERIALS & ASSIGNING SECTIONS ===")

for part_name in model.parts.keys():

    m = re.search(r"region_(\d+)", part_name)
    if not m:
        continue

    region_id = m.group(1)  # "008"
    print("\nProcessing part:", part_name, "| region_id:", region_id)

    if region_id not in material_database:
        print("WARNING: region_id not found in CSV:", region_id)
        continue

    props = material_database[region_id]
    material_name = "material_" + region_id
    section_name  = "Section_" + region_id
    set_name      = "Set_" + region_id

    # --- Create/replace Material ---
    if material_name in model.materials.keys():
        del model.materials[material_name]

    model.Material(name=material_name)
    mat = model.materials[material_name]

    # Thermal props (keep/remove as needed)
    mat.Density(table=((props["Density"],),))
    mat.Conductivity(table=((props["Conductivity"],),))
    mat.Expansion(table=((props["Expansion"],),))
    mat.SpecificHeat(law=CONSTANTPRESSURE, table=((props["SpecificHeat"],),))

    # --- UMAT constants (must be flat floats) ---
    umat_consts = props["UMAT"]
    mat.UserMaterial(mechanicalConstants=umat_consts)

    # SDVs
    if props.get("Depvar", None) is not None:
        mat.Depvar(n=int(props["Depvar"]))

    print("Created UMAT material:", material_name, "| #const =", len(umat_consts))

    # --- Create/replace Section ---
    if section_name in model.sections.keys():
        del model.sections[section_name]
    model.HomogeneousSolidSection(name=section_name, material=material_name, thickness=None)

    # --- Assign section to all cells ---
    p = model.parts[part_name]

    if len(p.cells) > 0:
        if set_name in p.sets.keys():
            del p.sets[set_name]
        region = p.Set(cells=p.cells, name=set_name)

        p.SectionAssignment(
            region=region,
            sectionName=section_name,
            offset=0.0,
            offsetType=MIDDLE_SURFACE,
            offsetField='',
            thicknessAssignment=FROM_SECTION
        )
        print("Assigned section:", section_name, "to:", part_name)
    else:
        print("WARNING: Part has no cells (surface-only import?):", part_name)

print("\n===== COMPLETE =====")

a = model.rootAssembly
a.DatumCsysByDefault(CARTESIAN)

# ---- USER CONTROLS (set as you need)
MERGED_PART_NAME     = "Part-1"      # merged part name
MERGED_INSTANCE_NAME = "Matrix-1"    # you want Matrix-1 == Part-1 instance
SUPPRESS_ORIGINALS   = True          # True = original region instances suppressed
mesh_size            = 2.5           # your size=2.5
deviationFactor      = 0.1
minSizeFactor        = 0.1

# ---- Collect region parts (region_###)
region_part_names = []
for pname in model.parts.keys():
    if re.match(r"^region_\d+$", pname):
        region_part_names.append(pname)

region_part_names.sort(key=lambda x: int(x.split("_")[1]))  # region_001,002,...

if len(region_part_names) == 0:
    raise RuntimeError("No parts found matching pattern: region_###")

# ---- Clean old instances of region_### if they exist
for iname in list(a.instances.keys()):
    if re.match(r"^region_\d+-\d+$", iname) or re.match(r"^region_\d+-1$", iname):
        try:
            del a.instances[iname]
        except:
            pass

# ---- Create instances in a loop
inst_names = []
for pname in region_part_names:
    iname = pname + "-1"
    if iname in a.instances.keys():
        del a.instances[iname]
    a.Instance(name=iname, part=model.parts[pname], dependent=ON)
    inst_names.append(iname)

print("Created instances:", inst_names)

# ---- Delete previous merged part/instance if present
if MERGED_PART_NAME in model.parts.keys():
    try:
        del model.parts[MERGED_PART_NAME]
    except:
        pass

if MERGED_PART_NAME + "-1" in a.instances.keys():
    try:
        del a.instances[MERGED_PART_NAME + "-1"]
    except:
        pass

if MERGED_INSTANCE_NAME in a.instances.keys():
    try:
        del a.instances[MERGED_INSTANCE_NAME]
    except:
        pass

# ---- Boolean merge
merge_instances = tuple(a.instances[nm] for nm in inst_names)

mergedFeat = a.InstanceFromBooleanMerge(
    name=MERGED_PART_NAME,
    instances=merge_instances,
    domain=GEOMETRY,
    keepIntersections=ON,
    originalInstances=SUPPRESS if SUPPRESS_ORIGINALS else KEEP
)
a.regenerate()

# Rename the created merged instance to Matrix-1 robustly
# Abaqus typically creates an instance called "<MERGED_PART_NAME>-1"
default_merged_instance = MERGED_PART_NAME + "-1"

if default_merged_instance in a.instances.keys():
    # rename the instance feature key (most robust)
    a.features.changeKey(fromName=default_merged_instance, toName=MERGED_INSTANCE_NAME)

a.regenerate()
print("Merged instance is:", MERGED_INSTANCE_NAME)

# ============================================================
# MESH THE MERGED PART (Part-level meshing; robust)
# ============================================================
pMerged = model.parts[MERGED_PART_NAME]

# Seed
pMerged.seedPart(size=mesh_size, deviationFactor=deviationFactor, minSizeFactor=minSizeFactor)

# Mesh controls (TET + FREE)  --- this matches what you were doing earlier
if len(pMerged.cells) == 0:
    raise RuntimeError("Merged part has no cells. Check your imported geometry (surface-only?).")

pMerged.setMeshControls(regions=pMerged.cells, elemShape=TET, technique=FREE)

# Element types (your earlier set: C3D20R/C3D15/C3D10). Keep if you want 2nd-order tets.
elemType1 = mesh.ElemType(elemCode=C3D20R, elemLibrary=STANDARD)
elemType2 = mesh.ElemType(elemCode=C3D15,  elemLibrary=STANDARD)
elemType3 = mesh.ElemType(elemCode=C3D10,  elemLibrary=STANDARD)
pMerged.setElementType(regions=(pMerged.cells,), elemTypes=(elemType1, elemType2, elemType3))

# Generate mesh
pMerged.generateMesh()
a.regenerate()

print("Meshing complete for merged part:", MERGED_PART_NAME)

# ============================================================
# OPTIONAL: Create Step (if not exists)
# ============================================================
if "Loading" not in model.steps.keys():
    model.StaticStep(name="Loading", previous="Initial", initialInc=0.01, maxInc=0.1)

# ============================================================
# PBC GEOMETRY SETS (ONLY Matrix-1 instance)
# Creates:
#   - vertex_1..vertex_8 and vertex_#_nodes
#   - edge_1..edge_12 (geometry edge sets)
#   - face_Xmin, face_Xmax, face_Ymin, face_Ymax, face_Zmin, face_Zmax
# ============================================================

a1 = model.rootAssembly
vertexlist,edgelist,matrix_facelist=[],[],[]
matrixInstance=a1.instances['Matrix-1']
vertices,edges=matrixInstance.vertices,matrixInstance.edges
vertexlist.append(vertices)
edgelist.append(edges)
matrix_facelist.append(a1.instances['Matrix-1'].faces)

# ---- You MUST have these defined somewhere:
wid = 100
hei = 100
thk = 100

# Tolerances (scale with size)
tolX = 1e-3 * wid
tolY = 1e-3 * hei
tolZ = 1e-3 * thk if thk != 0 else 1e-6

# Detect and create vertex sets ( Geom sets & node sets)

a1 = model.rootAssembly
vlist=[]
vset=[]
vertexloc=[[-wid/2,-hei/2,-thk/2],[wid/2,-hei/2,-thk/2],[wid/2,-hei/2,thk/2],[-wid/2,-hei/2,thk/2],[-wid/2,hei/2,-thk/2],[wid/2,hei/2,thk/2],[wid/2,hei/2,thk/2],[-wid/2,hei/2,thk/2]]
for i in range(len(vertexloc)):	
	for j in range(len(vertexlist)):
		flag = vertexlist[j].findAt(((vertexloc[i][0],vertexloc[i][1],vertexloc[i][2]),))	
		if flag:
			vlist.append(flag)
			setname='vertex_'+str(i+1)
			a1.Set(vertices=flag,name=setname)  
			set = a1.Set(name=setname+str('_nodes'), nodes=a1.sets[setname].nodes)
			vset.append(set)
            
# Create a list of edges
# edge.pointOn has been looped over all edges because edge.getByBoundingBox works inconsistantly with edges

edge1,edge2,edge3,edge4,edge5,edge6,edge7,edge8,edge9,edge10,edge11,edge12=[],[],[],[],[],[],[],[],[],[],[],[]
for i in range(len(edgelist)):
	curedges=edgelist[i]
	for j in range(len(curedges)):
		curedge=curedges[j]
		point=curedge.pointOn[0]
		if ((abs(point[0] -(-wid/2)) < 0.001*wid ) and (abs(point[1] - (-hei/2)) < 0.001*hei)):		#Edge 1
			edge1=edge1+[curedges.findAt(((point[0],point[1],point[2]),),)]
		if ((abs(point[0] - wid/2) < 0.001*wid) and (abs(point[1] - (-hei/2)) < 0.001*hei)):		#Edge 2
			edge2=edge2+[curedges.findAt(((point[0],point[1],point[2]),),)]
		if ((abs(point[0] - wid/2) < 0.001*wid) and (abs(point[1] - hei/2) < 0.001*hei)):		#Edge 3
			edge3=edge3+[curedges.findAt(((point[0],point[1],point[2]),),)]
		if ((abs(point[0] -( -wid/2)) < 0.001*wid) and (abs(point[1] - hei/2) < 0.001*hei)):		#Edge 4
			edge4=edge4+[curedges.findAt(((point[0],point[1],point[2]),),)]
		#
		if ((abs(point[0] -( -wid/2)) < 0.001*wid) and (abs(point[2] - (-thk/2)) <0.001*thk)):			#Edge 5
			edge5=edge5+[curedges.findAt(((point[0],point[1],point[2]),),)]
		if ((abs(point[0] - wid/2) < 0.001*wid) and (abs(point[2] - (-thk/2)) <0.001*thk)):			#Edge 6
			edge6=edge6+[curedges.findAt(((point[0],point[1],point[2]),),)]
		if ((abs(point[0] - wid/2) < 0.001*wid) and (abs(point[2] - (thk/2))<0.001*thk)):			#Edge 7
			edge7=edge7+[curedges.findAt(((point[0],point[1],point[2]),),)]
		if ((abs(point[0] -( -wid/2)) < 0.001*wid) and (abs(point[2] - (thk/2)) < 0.001*thk)):			#Edge 8
			edge8=edge8+[curedges.findAt(((point[0],point[1],point[2]),),)]
		#
		if ((abs(point[1] -( -hei/2)) < 0.001*hei) and (abs(point[2] - (-thk/2)) <0.001*thk)):				#Edge 9
			edge9=edge9+[curedges.findAt(((point[0],point[1],point[2]),),)]
		if ((abs(point[1] - hei/2) < 0.001*hei) and (abs(point[2] - (-thk/2)) <0.001*thk)):				#Edge 10
			edge10=edge10+[curedges.findAt(((point[0],point[1],point[2]),),)]
		if ((abs(point[1] - hei/2) < 0.001*hei) and (abs(point[2] - (thk/2)) <0.001*thk )):				#Edge 11
			edge11=edge11+[curedges.findAt(((point[0],point[1],point[2]),),)]
		if ((abs(point[1] - (-hei/2)) < 0.001*hei) and (abs(point[2] - (thk/2)) <0.001*thk)):				#Edge 12
			edge12=edge12+[curedges.findAt(((point[0],point[1],point[2]),),)]


a1.Set(edges=edge1, name='edge_1')
a1.Set(nodes=a1.sets['edge_1'].nodes, name = 'edge_1_nodes')
a1.Set(edges=edge2, name='edge_2')
a1.Set(nodes=a1.sets['edge_2'].nodes, name = 'edge_2_nodes')
a1.Set(edges=edge3, name='edge_3')
a1.Set(nodes=a1.sets['edge_3'].nodes, name = 'edge_3_nodes')
a1.Set(edges=edge4, name='edge_4')
a1.Set(nodes=a1.sets['edge_4'].nodes, name = 'edge_4_nodes')
a1.Set(edges=edge5, name='edge_5')
a1.Set(nodes=a1.sets['edge_5'].nodes, name = 'edge_5_nodes')
a1.Set(edges=edge6, name='edge_6')
a1.Set(nodes=a1.sets['edge_6'].nodes, name = 'edge_6_nodes')
a1.Set(edges=edge7, name='edge_7')
a1.Set(nodes=a1.sets['edge_7'].nodes, name = 'edge_7_nodes')
a1.Set(edges=edge8, name='edge_8')
a1.Set(nodes=a1.sets['edge_8'].nodes, name = 'edge_8_nodes')
a1.Set(edges=edge9, name='edge_9')
a1.Set(nodes=a1.sets['edge_9'].nodes, name = 'edge_9_nodes')
a1.Set(edges=edge10, name='edge_10')
a1.Set(nodes=a1.sets['edge_10'].nodes, name = 'edge_10_nodes')
a1.Set(edges=edge11, name='edge_11')
a1.Set(nodes=a1.sets['edge_11'].nodes, name = 'edge_11_nodes')
a1.Set(edges=edge12, name='edge_12')
a1.Set(nodes=a1.sets['edge_12'].nodes, name = 'edge_12_nodes')

# Detect and create face sets ( Geom sets & node sets) for Matrix

XinfF=XsupF=YinfF=YsupF=ZinfF=ZsupF=[]
for i in range(len(matrix_facelist)):
	cur_faces=matrix_facelist[i]
	XinfF=XinfF+[cur_faces.getByBoundingBox(-0.51*wid,-hei,-thk,-0.49*wid,hei,thk)]	#Xinf
	XsupF=XsupF+[cur_faces.getByBoundingBox(0.49*wid,-hei,-thk,0.51*wid,hei,thk)]	#Xsup
	YinfF=YinfF+[cur_faces.getByBoundingBox(-wid,-0.51*hei,-thk,wid,-0.49*hei,thk)] 	#Yinf
	YsupF=YsupF+[cur_faces.getByBoundingBox(-wid,0.49*hei,-thk,wid,0.51*hei,thk)] 	#Ysup
	ZinfF=ZinfF+[cur_faces.getByBoundingBox(-wid,-hei,-0.51*thk,wid,hei,-0.49*thk)] 	#Zinf
	ZsupF=ZsupF+[cur_faces.getByBoundingBox(-wid,-hei,0.49*thk,wid,hei,0.51*thk)]	#Zsup

a1.Set(faces=XinfF, name='X_inf_face')
a1.Set(nodes=a1.sets['X_inf_face'].nodes, name = 'X_inf_face_nodes')
a1.Set(faces=XsupF, name='X_sup_face')
a1.Set(nodes=a1.sets['X_sup_face'].nodes, name = 'X_sup_face_nodes')
a1.Set(faces=YinfF, name='Y_inf_face')
a1.Set(nodes=a1.sets['Y_inf_face'].nodes, name = 'Y_inf_face_nodes')
a1.Set(faces=YsupF, name='Y_sup_face')
a1.Set(nodes=a1.sets['Y_sup_face'].nodes, name = 'Y_sup_face_nodes')
a1.Set(faces=ZinfF, name='Z_inf_face')
a1.Set(nodes=a1.sets['Z_inf_face'].nodes, name = 'Z_inf_face_nodes')
a1.Set(faces=ZsupF, name='Z_sup_face')
a1.Set(nodes=a1.sets['Z_sup_face'].nodes, name = 'Z_sup_face_nodes')


vertex_1 = a1.sets['vertex_1_nodes']
vertex_2 = a1.sets['vertex_2_nodes']
vertex_3 = a1.sets['vertex_3_nodes']
vertex_4 = a1.sets['vertex_4_nodes']
vertex_5 = a1.sets['vertex_5_nodes']
vertex_6 = a1.sets['vertex_6_nodes']
vertex_7 = a1.sets['vertex_7_nodes']
vertex_8 = a1.sets['vertex_8_nodes']

edge_1 = a1.sets['edge_1_nodes']
edge_2 = a1.sets['edge_2_nodes']
edge_3 = a1.sets['edge_3_nodes']
edge_4 = a1.sets['edge_4_nodes']
edge_5 = a1.sets['edge_5_nodes']
edge_6 = a1.sets['edge_6_nodes']
edge_7 = a1.sets['edge_7_nodes']
edge_8 = a1.sets['edge_8_nodes']
edge_9 = a1.sets['edge_9_nodes']
edge_10 = a1.sets['edge_10_nodes']
edge_11 = a1.sets['edge_11_nodes']
edge_12 = a1.sets['edge_12_nodes']

x_back  = a1.sets['X_inf_face_nodes']
x_front = a1.sets['X_sup_face_nodes']
y_back  = a1.sets['Y_inf_face_nodes']
y_front = a1.sets['Y_sup_face_nodes']
z_back  = a1.sets['Z_inf_face_nodes']
z_front = a1.sets['Z_sup_face_nodes']

# ------------------------------------------------------------
# KUBC: 6 faces, 6 load cases (exx, eyy, ezz, gxy, gxz, gyz)
# ------------------------------------------------------------
# ---- Length scale (cube edge / RVE width) ----
L = wid  # keep your variable name
# ---- Output requests ----
model.fieldOutputRequests['F-Output-1'].setValues(variables=('S','U','IVOL','E'))
# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
def deleteBC_if_exists(bcname):
    if bcname in model.boundaryConditions.keys():
        del model.boundaryConditions[bcname]

def apply_sym_bc(bcname, region, normal):
    """
    normal = 'X' / 'Y' / 'Z'
    Applies symmetry in the normal direction only:
      XSYMM -> u1 = 0
      YSYMM -> u2 = 0
      ZSYMM -> u3 = 0
    """
    deleteBC_if_exists(bcname)
    u1 = UNSET; u2 = UNSET; u3 = UNSET
    if normal.upper() == 'X': u1 = 0.0
    if normal.upper() == 'Y': u2 = 0.0
    if normal.upper() == 'Z': u3 = 0.0

    model.DisplacementBC(
        name=bcname, createStepName='Loading', region=region,
        u1=u1, u2=u2, u3=u3, ur1=UNSET, ur2=UNSET, ur3=UNSET,
        amplitude=UNSET, fixed=OFF, distributionType=UNIFORM,
        fieldName='', localCsys=None
    )

def apply_disp_bc(bcname, region, dof, value):
    """
    dof = 1/2/3 corresponding to u1/u2/u3
    Applies only that dof, leaves others UNSET.
    """
    deleteBC_if_exists(bcname)
    u1 = UNSET; u2 = UNSET; u3 = UNSET
    if dof == 1: u1 = value
    if dof == 2: u2 = value
    if dof == 3: u3 = value

    model.DisplacementBC(
        name=bcname, createStepName='Loading', region=region,
        u1=u1, u2=u2, u3=u3, ur1=UNSET, ur2=UNSET, ur3=UNSET,
        amplitude=UNSET, fixed=OFF, distributionType=UNIFORM,
        fieldName='', localCsys=None
    )

def apply_common_symmetry_except_loaded(loaded_case):
    """
    Put symmetry BCs on faces as needed.
    For tensile in X: x_back XSYM, y_back&front YSYM, z_back&front ZSYM
    For tensile in Y: y_back YSYM, x_back&front XSYM, z_back&front ZSYM
    For tensile in Z: z_back ZSYM, x_back&front XSYM, y_back&front YSYM
    For shear:
      g12: x_back XSYM, y_back&front YSYM, z_back&front ZSYM (then apply u1 on y_front)
      g13: x_back XSYM, y_back&front YSYM, z_back&front ZSYM (then apply u1 on z_front)
      g23: x_back&front XSYM, y_back YSYM, z_back&front ZSYM (then apply u2 on z_front)
    """
    # First delete any old symmetry BCs (optional cleanup)
    for nm in ['BC_XBACK_XSYM','BC_XFRONT_XSYM','BC_YBACK_YSYM','BC_YFRONT_YSYM','BC_ZBACK_ZSYM','BC_ZFRONT_ZSYM']:
        deleteBC_if_exists(nm)

    if loaded_case == 'EXX':
        apply_sym_bc('BC_XBACK_XSYM',  x_back,  'X')
        apply_sym_bc('BC_YBACK_YSYM',  y_back,  'Y')
        apply_sym_bc('BC_YFRONT_YSYM', y_front, 'Y')
        apply_sym_bc('BC_ZBACK_ZSYM',  z_back,  'Z')
        apply_sym_bc('BC_ZFRONT_ZSYM', z_front, 'Z')

    elif loaded_case == 'EYY':
        apply_sym_bc('BC_YBACK_YSYM',  y_back,  'Y')
        apply_sym_bc('BC_XBACK_XSYM',  x_back,  'X')
        apply_sym_bc('BC_XFRONT_XSYM', x_front, 'X')
        apply_sym_bc('BC_ZBACK_ZSYM',  z_back,  'Z')
        apply_sym_bc('BC_ZFRONT_ZSYM', z_front, 'Z')

    elif loaded_case == 'EZZ':
        apply_sym_bc('BC_ZBACK_ZSYM',  z_back,  'Z')
        apply_sym_bc('BC_XBACK_XSYM',  x_back,  'X')
        apply_sym_bc('BC_XFRONT_XSYM', x_front, 'X')
        apply_sym_bc('BC_YBACK_YSYM',  y_back,  'Y')
        apply_sym_bc('BC_YFRONT_YSYM', y_front, 'Y')

    elif loaded_case == 'G12':   # xy shear
        apply_sym_bc('BC_XBACK_XSYM',  x_back,  'X')
        apply_sym_bc('BC_YBACK_YSYM',  y_back,  'Y')
        apply_sym_bc('BC_YFRONT_YSYM', y_front, 'Y')   # keeps u2=0, but u1 is free -> OK
        apply_sym_bc('BC_ZBACK_ZSYM',  z_back,  'Z')
        apply_sym_bc('BC_ZFRONT_ZSYM', z_front, 'Z')

    elif loaded_case == 'G13':   # xz shear
        apply_sym_bc('BC_XBACK_XSYM',  x_back,  'X')
        apply_sym_bc('BC_YBACK_YSYM',  y_back,  'Y')
        apply_sym_bc('BC_YFRONT_YSYM', y_front, 'Y')
        apply_sym_bc('BC_ZBACK_ZSYM',  z_back,  'Z')
        apply_sym_bc('BC_ZFRONT_ZSYM', z_front, 'Z')   # keeps u3=0, but u1 is free -> OK

    elif loaded_case == 'G23':   # yz shear
        apply_sym_bc('BC_ZBACK_ZSYM',  z_back,  'Z')
        apply_sym_bc('BC_ZFRONT_ZSYM', z_front, 'Z')   # keeps u3=0, but u2 is free -> OK
        apply_sym_bc('BC_YBACK_YSYM',  y_back,  'Y')    # keeps u2=0 on y_back -> reference
        apply_sym_bc('BC_XBACK_XSYM',  x_back,  'X')
        apply_sym_bc('BC_XFRONT_XSYM', x_front, 'X')


# ------------------------------------------------------------
# Load cases: 6 jobs
#   EXX: apply u1 on x_front
#   EYY: apply u2 on y_front
#   EZZ: apply u3 on z_front
#   G12: apply u1 on y_front
#   G13: apply u1 on z_front
#   G23: apply u2 on z_front
# ------------------------------------------------------------
cases = [
    ('EXX', 0.01),  # e11
    ('EYY', 0.01),  # e22
    ('EZZ', 0.01),  # e33
    ('G12', 0.01),  # g12
    ('G13', 0.01),  # g13
    ('G23', 0.01),  # g23
]

for i, (caseName, eps) in enumerate(cases):

    # Apply symmetry set for this case
    apply_common_symmetry_except_loaded(caseName)

    # Remove any previous "load" BC
    deleteBC_if_exists('BC_LOAD')

    # Apply the actual loading displacement
    disp = eps * L

    if caseName == 'EXX':
        # Your requested pattern:
        # x_back XSYMM, y faces YSYMM, z faces ZSYMM, x_front u1 = disp
        apply_disp_bc('BC_LOAD', x_front, 1, disp)

    elif caseName == 'EYY':
        apply_disp_bc('BC_LOAD', y_front, 2, disp)

    elif caseName == 'EZZ':
        apply_disp_bc('BC_LOAD', z_front, 3, disp)

    elif caseName == 'G12':
        # xy shear: u1 on y_front
        apply_disp_bc('BC_LOAD', y_front, 1, disp)

    elif caseName == 'G13':
        # xz shear: u1 on z_front
        apply_disp_bc('BC_LOAD', z_front, 1, disp)

    elif caseName == 'G23':
        # yz shear: u2 on z_front
        apply_disp_bc('BC_LOAD', z_front, 2, disp)

    # Create job
    jobname = 'Mesoscale_' + caseName
    mdb.Job(
        name=jobname, model='Model-1', description='KUBC case: '+caseName,
        type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0,
        queue=None, memory=50, memoryUnits=PERCENTAGE,
        getMemoryFromAnalysis=False, explicitPrecision=SINGLE,
        nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF,
        contactPrint=OFF, historyPrint=OFF, userSubroutine='',
        scratch='', multiprocessingMode=DEFAULT, numCpus=5, numDomains=5
    )
    mdb.jobs[jobname].writeInput(consistencyChecking=OFF)