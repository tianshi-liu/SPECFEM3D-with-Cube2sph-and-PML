import os
import numpy as np
R_EARTH = 6371000.0
METER_PER_DEG = R_EARTH * np.pi / 180.0
x_size_deg = 22.0
y_size_deg = 22.0
truncation_depth = 500.0 # in km
nelem_x = 160
nelem_y = 160
nelem_z = [4,         7,         5,        3,        1,        1] # number of elements in each vertical region
top_z   = [-400000.0, -220000.0, -80000.0, -24400.0, -15000.0, 0.0] # z coordinate at the top of each vertical region
fn_int  = ['interface_400.dat', 'interface_220.dat', 'interface_80.dat', 'interface_moho.dat', 'interface_15.dat', 'interface_topo.dat']
idoubling=[                                 3] # doubling at the top of this vertical region
ilayer   =[                                 3,                            6] # layering at the top of this vertical region
#           rho     Vp        Vs        Q_Kappa Q_mu  ani_flag dom_id 
prop     =[[3849.8, 9645.870, 5224.210, 9999.0, 40.0, 0,       2],
           [2600.0, 5800.000, 3200.000, 9999.0, 50.0, 0,       2]]
nelem_cpml_x = 4
nelem_cpml_y = 4
nelem_cpml_z = 2

mesh_dir = 'DATA/meshfem3D_files'
fn_mesh = 'Mesh_Par_file'
fn_interface = 'interfaces.dat'
fn = os.path.join(mesh_dir, fn_mesh)
f = open(fn, 'w')

f.write("#-----------------------------------------------------------\n")
f.write("#\n")
f.write("# Meshing input parameters\n")
f.write("#\n")
f.write("#-----------------------------------------------------------\n\n")
f.write("# coordinates of mesh block in latitude/longitude and depth in km\n")
f.write(f"LATITUDE_MIN                    = -{(METER_PER_DEG * y_size_deg / 2.0):.4f}\n")
f.write(f"LATITUDE_MAX                    = {(METER_PER_DEG * y_size_deg / 2.0):.4f}\n")
f.write(f"LONGITUDE_MIN                    = -{(METER_PER_DEG * x_size_deg / 2.0):.4f}\n")
f.write(f"LONGITUDE_MAX                    = {(METER_PER_DEG * x_size_deg / 2.0):.4f}\n")
f.write(f"DEPTH_BLOCK_KM                  = {truncation_depth:.2f}\n")
f.write("UTM_PROJECTION_ZONE             = 11\n") # will be ignored
f.write("SUPPRESS_UTM_PROJECTION         = .true.\n\n")

f.write("# file that contains the interfaces of the model / mesh\n")
f.write(f"INTERFACES_FILE                 = {fn_interface}\n\n")

f.write("# file that contains the cavity\n")
f.write("CAVITY_FILE                     = no_cavity.dat\n\n") # not supported, will be ignored

f.write("# number of elements at the surface along edges of the mesh at the surface\n")
f.write("# (must be 8 * multiple of NPROC below if mesh is not regular and contains mesh doublings)\n")
f.write("# (must be multiple of NPROC below if mesh is regular)\n")
f.write(f"NEX_XI                          = {nelem_x}\n")
f.write(f"NEX_ETA                         = {nelem_y}\n\n")

f.write("# number of MPI processors along xi and eta (can be different)\n")
f.write("NPROC_XI                        = 1\n")
f.write("NPROC_ETA                       = 1\n\n") # will be ignored

f.write("#-----------------------------------------------------------\n")
f.write("#\n")
f.write("# Doubling layers\n")
f.write("#\n")
f.write("#-----------------------------------------------------------\n\n")

f.write("# Regular/irregular mesh\n")
f.write("USE_REGULAR_MESH                = .false.\n")
f.write("# Only for irregular meshes, number of doubling layers and their position\n")
f.write(f"NDOUBLINGS                      = {len(idoubling)}\n")
f.write("# NZ_DOUBLING_1 is the parameter to set up if there is only one doubling layer\n")
f.write("# (more doubling entries can be added if needed to match NDOUBLINGS value)\n")
for i in range(0, len(idoubling)):
  f.write(f"NZ_DOUBLING_{i+1}                   = {np.sum(np.array(nelem_z, dtype=int)[0:idoubling[i]+1]):d}\n")
f.write(f"NZ_DOUBLING_{len(idoubling)+1}                   = 0\n\n")

f.write("#-----------------------------------------------------------\n")
f.write("#\n")
f.write("# Visualization\n")
f.write("#\n")
f.write("#-----------------------------------------------------------\n\n")

f.write("# create mesh files for visualisation or further checking\n")
f.write("CREATE_ABAQUS_FILES             = .false.\n")
f.write("CREATE_DX_FILES                 = .false.\n")
f.write("CREATE_VTK_FILES                = .true.\n\n")

f.write("# path to store the databases files\n")
f.write("LOCAL_PATH                      = ./DATABASES_MPI\n\n")

f.write("#-----------------------------------------------------------\n")
f.write("#\n")
f.write("# CPML\n")
f.write("#\n")
f.write("#-----------------------------------------------------------\n\n")

f.write("# CPML perfectly matched absorbing layers, in m\n")
f.write(f"THICKNESS_OF_X_PML              = {(x_size_deg * METER_PER_DEG / nelem_x * nelem_cpml_x):.4f}\n")
f.write(f"THICKNESS_OF_Y_PML              = {(y_size_deg * METER_PER_DEG / nelem_y * nelem_cpml_y):.4f}\n")
f.write(f"THICKNESS_OF_Z_PML              = {((top_z[0] + truncation_depth * 1000.0) / nelem_z[0] * nelem_cpml_z):.4f}\n")

f.write("#-----------------------------------------------------------\n")
f.write("#\n")
f.write("# Domain materials\n")
f.write("#\n")
f.write("#-----------------------------------------------------------\n\n")

f.write("# number of materials\n")
f.write(f"NMATERIALS                      = {len(ilayer)}\n")
f.write("# define the different materials in the model as:\n")
f.write("# #material_id  #rho  #vp  #vs  #Q_Kappa  #Q_mu  #anisotropy_flag  #domain_id\n")
f.write("#     Q_Kappa          : Q_Kappa attenuation quality factor\n")
f.write("#     Q_mu             : Q_mu attenuation quality factor\n")
f.write("#     anisotropy_flag  : 0 = no anisotropy / 1,2,... check the implementation in file aniso_model.f90\n")
f.write("#     domain_id        : 1 = acoustic / 2 = elastic\n")
for i in range(0, len(ilayer)):
  f.write(f"{i+1} {prop[i][0]:.1f} {prop[i][1]:.3f} {prop[i][2]:.3f} {prop[i][3]:.0f} {prop[i][4]:.1f} {prop[i][5]} {prop[i][6]}\n")
f.write("\n")

f.write("#-----------------------------------------------------------\n")
f.write("#\n")
f.write("# Domain regions\n")
f.write("#\n")
f.write("#-----------------------------------------------------------\n\n")

f.write("# number of regions\n")
f.write(f"NREGIONS                        = {len(ilayer)}\n")
f.write("# define the different regions of the model as :\n")
f.write("#NEX_XI_BEGIN  #NEX_XI_END  #NEX_ETA_BEGIN  #NEX_ETA_END  #NZ_BEGIN #NZ_END  #material_id\n")
f.write(f"1    {nelem_x}    1    {nelem_y}    1    {np.sum(np.array(nelem_z)[0:ilayer[0]+1]):d} 1\n")
for i in range(0, len(ilayer)-1):
  f.write(f"1    {nelem_x}    1    {nelem_y}    {np.sum(np.array(nelem_z)[0:ilayer[i]+1])+1:d}    {np.sum(np.array(nelem_z)[0:ilayer[i+1]+1]):d} {i+2}\n")

f.close()

fn = os.path.join(mesh_dir, fn_interface)
f = open(fn, 'w')
f.write("# number of interfaces\n")
f.write(f" {len(nelem_z)}\n")
f.write("#\n")
f.write("# We describe each interface below, structured as a 2D-grid, with several parameters : \n")
f.write("# number of points along XI and ETA, minimal XI ETA coordinates\n")
f.write("# and spacing between points which must be constant.\n")
f.write("# Then the records contain the Z coordinates of the NXI x NETA points.\n")
f.write("#\n")
f.write("# SUPPRESS_UTM_PROJECTION  NXI  NETA LONG_MIN   LAT_MIN    SPACING_XI SPACING_ETA\n")
for i in range(0, len(fn_int)):
  f.write(f" .true.  2  2  -{(METER_PER_DEG * x_size_deg / 2.0):.4f}  -{(METER_PER_DEG * y_size_deg / 2.0):.4f}  {(METER_PER_DEG * x_size_deg ):.4f}  {(METER_PER_DEG * y_size_deg ):.4f}\n {fn_int[i]}\n")

f.write("# for each layer, we give the number of spectral elements in the vertical direction\n")
for i in range(0, len(nelem_z)):
  f.write(f"{nelem_z[i]}\n")
f.close()

for i in range(0, len(fn_int)):
  f = open(os.path.join(mesh_dir, fn_int[i]), 'w')
  f.write(f"{top_z[i]}\n"*4)
  f.close()
