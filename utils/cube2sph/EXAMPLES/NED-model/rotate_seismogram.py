"""
Python script to rotate seismograms between the Cartesian system and the Geographic system
Tianshi Liu - 2022.03.24

command line arguments:
fn_matrix: name of the rotation matrix file produced by the write_station_file program. Rows: (N, E, Z) in geographic system; Columns: (x, y, z) in Cartesian system
rotate: rotation scheme. "XYZ->NEZ" for Cartesian to geographic, "XYZ<-NEZ" for geographic back to Cartesian (this might be useful to transform adjoint sources). If some component is absent, use 0 for that component. For example, if I only use the vertical component, I can use "XYZ->00Z" and "XYZ<-00Z".
from_dir, to_dir: input and output directory
from_template, to_template: file name template for input and output. ${nt}, ${sta}, ${comp} should be used to stand for network, station and component. Note that these two arguments should use single quotation marks instead of double quotation marks in bash, otherwise bash will try to interpret ${}.

Examples:
(1)
> python3 rotate_seismogram.py --fn_matrix="rotation_nu" --rotate="XYZ->NEZ" --from_dir="OUTPUT_FILES" --to_dir="OUTPUT_FILES_sph" --from_template='${nt}.${sta}.BX${comp}.semd' --to_template='${nt}.${sta}.BX${comp}.sem.ascii'

inputs are in Cartesian (output by SPECFEM3D+Cube2sph), outputs are in geographic (output by e.g., SPECFEM3D_globe):
reading from OUTPUT_FILES/TS.TS11.BXX.semd
reading from OUTPUT_FILES/TS.TS11.BXY.semd
reading from OUTPUT_FILES/TS.TS11.BXZ.semd
writing to $OUTPUT_FILES_sph/TS.TS11.BXN.sem.ascii
writing to $OUTPUT_FILES_sph/TS.TS11.BXE.sem.ascii
writing to $OUTPUT_FILES_sph/TS.TS11.BXZ.sem.ascii
...

(2)
> python3 rotate_seismogram.py --fn_matrix="rotation_nu" --rotate="XYZ<-NEZ" --from_dir="OUTPUT_FILES_sph" --to_dir="OUTPUT_FILES_cart" --to_template='${nt}.${sta}.BX${comp}.semd' --from_template='${nt}.${sta}.BX${comp}.sem.ascii'

transform back:
reading from OUTPUT_FILES_sph/TS.TS11.BXN.sem.ascii
reading from OUTPUT_FILES_sph/TS.TS11.BXE.sem.ascii
reading from OUTPUT_FILES_sph/TS.TS11.BXZ.sem.ascii
writing to $OUTPUT_FILES_cart/TS.TS11.BXX.semd
writing to $OUTPUT_FILES_cart/TS.TS11.BXY.semd
writing to $OUTPUT_FILES_cart/TS.TS11.BXZ.semd
...

(3)
> python3 rotate_seismogram.py --fn_matrix="rotation_nu" --rotate="XYZ->00Z" --from_dir="OUTPUT_FILES" --to_dir="OUTPUT_FILES_sph" --from_template='${nt}.${sta}.BX${comp}.semd' --to_template='${nt}.${sta}.BX${comp}.sem.ascii'

ignoring N and E components
"""
import mpi4py.rc
mpi4py.rc.threads = False
import numpy as np
import argparse
import os
from string import Template
import subprocess
from mpi4py import MPI

parser = argparse.ArgumentParser()
parser.add_argument('--fn_matrix', dest='fn_matrix', type=str, metavar='FN_MATRIX')
parser.add_argument('--rotate', dest='rotate', type=str, metavar='ROTATE')
parser.add_argument('--from_dir', dest='from_dir', type=str, metavar='FROM_DIR')
parser.add_argument('--from_template', dest='from_template', type=str, metavar='FROM_TEMPLATE')
parser.add_argument('--to_dir', dest='to_dir', type=str, metavar='TO_DIR')
parser.add_argument('--to_template', dest='to_template', type=str, metavar='TO_TEMPLATE')
#parser.add_argument('--nsteps', dest='nsteps', type=int, metavar='NSTEPS')
args = parser.parse_args()

fn_matrix = args.fn_matrix
rotate = args.rotate
from_dir = args.from_dir
from_template = Template(args.from_template)
to_dir = args.to_dir
to_template = Template(args.to_template)

comp_left = rotate[0:3]
comp_right = rotate[5:8]
if (rotate[3:5] == '->'):
  forward = True
  from_comp, to_comp = comp_left, comp_right
elif (rotate[3:5] == '<-'):
  forward = False
  from_comp, to_comp = comp_right, comp_left
else:
  raise ValueError('invalid rotate parameter')

#arr = np.zeros(shape(nsteps, 2), dtype=float)
nu = np.zeros(shape=(3,3), dtype=float)
with open(fn_matrix, 'r') as rot_sta:
  lines = rot_sta.readlines()
#--MPI--
comm = MPI.COMM_WORLD
myrank = comm.Get_rank()
nproc = comm.Get_size()
nrec = len(lines)//4
nrec_local = nrec // nproc
if (myrank < (nrec % nproc)): nrec_local = nrec_local + 1
for i_sta_local in range(0, nrec_local):
  i_sta = i_sta_local * nproc + myrank
  line_segs = lines[i_sta*4].strip().split(' ')
  line_segs = [_ for _ in line_segs if _ != '']
  nt = line_segs[0]
  sta = line_segs[1]
  for i_comp in range(0, 3):
    line_segs = lines[i_sta*4+i_comp+1].strip().split(' ')
    line_segs = [_ for _ in line_segs if _ != '']
    nu[i_comp,0] = float(line_segs[0])
    nu[i_comp,1] = float(line_segs[1])
    nu[i_comp,2] = float(line_segs[2])
  missing_file = False
  nstep = -1
  for i_comp in range(0, 3):
    if (from_comp[i_comp] == '0'):
      continue 
    fn = os.path.join(from_dir, 
             from_template.substitute(nt=nt, sta=sta, comp=from_comp[i_comp]))
    if (not os.path.isfile(fn)):
      #print(f"{fn} does not exist but required by rotation, skipping this station\n")
      missing_file = True
      break
    if (nstep < 0):
      nstep = int(subprocess.check_output("wc -l " + fn, shell=True).split()[0])
  if (missing_file):
    continue
  seis = np.zeros(shape=(nstep, 3), dtype=float)
  for i_comp in range(0, 3):
    if (from_comp[i_comp] == '0'):
      continue
    fn = os.path.join(from_dir, 
             from_template.substitute(nt=nt, sta=sta, comp=from_comp[i_comp]))
    #print(f"reading from {fn}")
    arr = np.loadtxt(fn)
    seis[:,i_comp] = arr[:,1]

  if (forward):
    seis = np.matmul(seis, np.transpose(nu))
  else:
    seis = np.matmul(seis, nu)

  for i_comp in range(0, 3):
    if (to_comp[i_comp] == '0'):
      continue
    fn = os.path.join(to_dir, 
             to_template.substitute(nt=nt, sta=sta, comp=to_comp[i_comp]))
    arr[:,1] = seis[:,i_comp]
    #print(f"writing to ${fn}")
    np.savetxt(fname=fn, X=arr, fmt='%11.6f%19.7E')
         
      
