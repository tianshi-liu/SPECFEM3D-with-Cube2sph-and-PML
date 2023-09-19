import numpy as np
import os

#mesh8_fn = 'MESH-default/mesh_file'
home_dir = '.'
name_suffix_8 = ''
name_suffix_27 = '_27'
mesh8_fn = os.path.join(home_dir, 'mesh_file' + name_suffix_8)
mesh27_fn = os.path.join(home_dir, 'mesh_file' + name_suffix_27)

with open(mesh8_fn, 'r') as f_mesh8:
    n_element = f_mesh8.readline()

n_element = int(n_element)

mesh_nodes = np.loadtxt(fname=mesh8_fn, dtype=int, skiprows=1)

if not n_element == len(mesh_nodes):
    print('element number not consistent')
    exit(1)
print('there are ' + str(n_element) + ' elements.')

face_dict = {}
edge_dict = {}

face_to_node = []  # node index starting from 1
face_used_time = []
element_face = []  # face index starting from 0

edge_to_node = []
edge_used_time = []
element_edge = []

face_each_element = [
    [1, 2, 3, 4],
    [1, 5, 6, 2],
    [2, 3, 7, 6],
    [4, 8, 7, 3],
    [1, 4, 8, 5],
    [5, 6, 7, 8]
]
edge_each_element = [
    [1, 2], [2, 3], [4, 3], [1, 4],
    [1, 5], [2, 6], [3, 7], [4, 8],
    [5, 6], [6, 7], [8, 7], [5, 8]
]
n_faces = 0
n_edges = 0
face_per_element = 6
edge_per_element = 12
max_face_element = 2
max_edge_element = 4

node_perm_face = [
    [0, 1, 2, 3],
    [1, 2, 3, 0],
    [2, 3, 0, 1],
    [3, 0, 1, 2],
    [3, 2, 1, 0],
    [0, 3, 2, 1],
    [1, 0, 3, 2],
    [2, 1, 0, 3]
]
node_perm_edge = [
    [0, 1],
    [1, 0]
]
print_every = 10000
print('counting faces...')
for line in mesh_nodes:
    i_element = line[0]
    _, r = divmod(i_element, print_every)
    if r == 0:
        print('element ' + str(i_element) + '\n')
    node = line
    this_element_face = []
    for each_face_of_element in face_each_element:
        this_face_node = node[each_face_of_element]
        this_face_id = None
        for perm in node_perm_face:
            this_face_id = face_dict.get(tuple(this_face_node[perm]), None)
            if this_face_id is not None:
                right_node_order = tuple(this_face_node[perm])
                break
        if this_face_id is None:
            this_face_id = n_faces
            n_faces = n_faces + 1
            face_dict[tuple(this_face_node)] = this_face_id
            face_used_time.append(1)
            face_to_node.append(this_face_node)
        else:
            face_used_time[this_face_id] += 1
            if face_used_time[this_face_id] == max_face_element:
                face_dict.pop(right_node_order)
        this_element_face.append(this_face_id)
    element_face.append(this_element_face)
print('there are ' + str(n_faces) + ' faces.')
len_dict = len(face_dict)
print('counting edges...')
for line in mesh_nodes:
    i_element = line[0]
    _, r = divmod(i_element, print_every)
    if r == 0:
        print('element ' + str(i_element) + '\n')
    node = line
    this_element_edge = []
    for each_edge_of_element in edge_each_element:
        this_edge_node = node[each_edge_of_element]
        this_edge_id = None
        for perm in node_perm_edge:
            this_edge_id = edge_dict.get(tuple(this_edge_node[perm]), None)
            if this_edge_id is not None:
                right_node_order = tuple(this_edge_node[perm])
                break
        if this_edge_id is None:
            this_edge_id = n_edges
            n_edges = n_edges + 1
            edge_dict[tuple(this_edge_node)] = this_edge_id
            edge_used_time.append(1)
            edge_to_node.append(this_edge_node)
        else:
            edge_used_time[this_edge_id] += 1
            if edge_used_time[this_edge_id] == max_edge_element:
                edge_dict.pop(right_node_order)
        this_element_edge.append(this_edge_id)
    element_edge.append(this_element_edge)
print('there are ' + str(n_edges) + ' edges.')
print('start writing files of HEX27 elements...')
mesh8_coords_fn = os.path.join(home_dir ,'nodes_coords_file' + name_suffix_8)
mesh27_coords_fn = os.path.join(home_dir, 'nodes_coords_file' + name_suffix_27)
with open(mesh8_coords_fn, 'r') as f_coord8:
    n_vertex = f_coord8.readline()
n_node_per_edge = 2
n_node_per_face = 4
n_node_per_element = 8
n_vertex = int(n_vertex)
vertex_coords = np.genfromtxt(fname=mesh8_coords_fn, dtype=[('inode', 'i4'), ('x', 'f8'), ('y', 'f8'), ('z', 'f8')], skip_header=1)
print('computing positions of edge centers...')
edge_center_coords = np.zeros(shape=(n_edges, ), dtype=[('inode', 'i4'), ('x', 'f8'), ('y', 'f8'), ('z', 'f8')])
edge_to_node = np.array(edge_to_node, dtype=int)
#x_temp = vertex_coords[:, 1]
#x_temp = x_temp[edge_to_node]
edge_center_coords['x'] = np.sum(vertex_coords['x'][edge_to_node - 1], axis=1) / n_node_per_edge
edge_center_coords['y'] = np.sum(vertex_coords['y'][edge_to_node - 1], axis=1) / n_node_per_edge
edge_center_coords['z'] = np.sum(vertex_coords['z'][edge_to_node - 1], axis=1) / n_node_per_edge
edge_center_offset = n_vertex + 1
edge_center_coords['inode'] = np.arange(start=edge_center_offset, stop=edge_center_offset + n_edges, step=1, dtype=int)


print('computing positions of face centers...')
face_center_coords = np.zeros(shape=(n_faces, ), dtype=[('inode', 'i4'), ('x', 'f8'), ('y', 'f8'), ('z', 'f8')])
face_to_node = np.array(face_to_node, dtype=int)
face_center_coords['x'] = np.sum(vertex_coords['x'][face_to_node - 1], axis=1) / n_node_per_face
face_center_coords['y'] = np.sum(vertex_coords['y'][face_to_node - 1], axis=1) / n_node_per_face
face_center_coords['z'] = np.sum(vertex_coords['z'][face_to_node - 1], axis=1) / n_node_per_face
face_center_offset = n_vertex + n_edges + 1
face_center_coords['inode'] = np.arange(start=face_center_offset, stop=face_center_offset + n_faces, step=1, dtype=int)

print('computing positions of element centers...')
element_center_coords = np.zeros(shape=(n_element, ), dtype=[('inode', 'i4'), ('x', 'f8'), ('y', 'f8'), ('z', 'f8')])
#face_to_node = np.array(face_to_node, dtype=int)
element_center_coords['x'] = np.sum(vertex_coords['x'][mesh_nodes[:, 1:] - 1], axis=1) / n_node_per_element
element_center_coords['y'] = np.sum(vertex_coords['y'][mesh_nodes[:, 1:] - 1], axis=1) / n_node_per_element
element_center_coords['z'] = np.sum(vertex_coords['z'][mesh_nodes[:, 1:] - 1], axis=1) / n_node_per_element
element_center_offset = n_vertex + n_edges + n_faces + 1
element_center_coords['inode'] = np.arange(start=element_center_offset, stop=element_center_offset + n_element, step=1, dtype=int)

n_nodes_27 = n_vertex + n_edges + n_faces + n_element
nodes_coords_27 = np.zeros(shape=(n_nodes_27, ), dtype=[('inode', 'i4'), ('x', 'f8'), ('y', 'f8'), ('z', 'f8')])
nodes_coords_27['inode'] = np.concatenate((vertex_coords['inode'], edge_center_coords['inode'], face_center_coords['inode'], element_center_coords['inode']), axis=0)
nodes_coords_27['x'] = np.concatenate((vertex_coords['x'], edge_center_coords['x'], face_center_coords['x'], element_center_coords['x']), axis=0)
nodes_coords_27['y'] = np.concatenate((vertex_coords['y'], edge_center_coords['y'], face_center_coords['y'], element_center_coords['y']), axis=0)
nodes_coords_27['z'] = np.concatenate((vertex_coords['z'], edge_center_coords['z'], face_center_coords['z'], element_center_coords['z']), axis=0)
np.savetxt(fname=mesh27_coords_fn, X=nodes_coords_27, header=str(n_nodes_27), comments='',
           fmt=['%12i', '%20.6f', '%20.6f', '%20.6f'])


print('constructing HEX27 elements...')
element_edge = np.array(element_edge, dtype=int)
element_face = np.array(element_face, dtype=int)
mesh27_nodes = np.concatenate((mesh_nodes, element_edge + edge_center_offset, element_face + face_center_offset, np.reshape(element_center_coords['inode'], newshape=(n_element, 1))), axis=1)
np.savetxt(fname=mesh27_fn, X=mesh27_nodes, header=str(n_element), comments='',
           fmt=['%-12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i',
                '%12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i',
                '%12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i'])

idx_bottom = [4, 3, 2, 1, 11, 10, 9, 12, 21]
idx_xmin = [1, 5, 8, 4, 13, 20, 16, 12, 25]
idx_xmax = [2, 3, 7, 6, 10, 15, 18, 14, 23]
idx_ymin = [1, 2, 6, 5, 9, 14, 17, 13, 22]
idx_ymax = [4, 8, 7, 3, 16, 19, 15, 11, 24]
idx_top = [5, 6, 7, 8, 17, 18, 19, 20, 26]
n_abs_total = 0

print('writing absorbing boundary files...')
abs_fn = os.path.join(home_dir, 'absorbing_surface_file_bottom' + name_suffix_8)
with open(abs_fn, 'r') as f_abs:
    n_abs = f_abs.readline()
n_abs = int(n_abs)
n_abs_total += n_abs
abs_nodes = np.loadtxt(fname=abs_fn, dtype=int, skiprows=1)
element_abs = abs_nodes[:, 0]
abs_nodes = np.concatenate((np.reshape(element_abs, newshape=(n_abs, 1)), mesh27_nodes[np.ix_(element_abs - 1, idx_bottom)]), axis=1)
abs_fn = os.path.join(home_dir, 'absorbing_surface_file_bottom' + name_suffix_27)
np.savetxt(fname=abs_fn, X=abs_nodes, header=str(n_abs), comments='',
           fmt=['%-12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i'])


abs_fn = os.path.join(home_dir, 'absorbing_surface_file_xmax' + name_suffix_8)
with open(abs_fn, 'r') as f_abs:
    n_abs = f_abs.readline()
n_abs = int(n_abs)
n_abs_total += n_abs
abs_nodes = np.loadtxt(fname=abs_fn, dtype=int, skiprows=1)
element_abs = abs_nodes[:, 0]
abs_nodes = np.concatenate((np.reshape(element_abs, newshape=(n_abs, 1)), mesh27_nodes[np.ix_(element_abs - 1, idx_xmax)]), axis=1)
abs_fn = os.path.join(home_dir, 'absorbing_surface_file_xmax' + name_suffix_27)
np.savetxt(fname=abs_fn, X=abs_nodes, header=str(n_abs), comments='',
           fmt=['%-12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i'])

abs_fn = os.path.join(home_dir, 'absorbing_surface_file_xmin' + name_suffix_8)
with open(abs_fn, 'r') as f_abs:
    n_abs = f_abs.readline()
n_abs = int(n_abs)
n_abs_total += n_abs
abs_nodes = np.loadtxt(fname=abs_fn, dtype=int, skiprows=1)
element_abs = abs_nodes[:, 0]
abs_nodes = np.concatenate((np.reshape(element_abs, newshape=(n_abs, 1)), mesh27_nodes[np.ix_(element_abs - 1, idx_xmin)]), axis=1)
abs_fn = os.path.join(home_dir, 'absorbing_surface_file_xmin' + name_suffix_27)
np.savetxt(fname=abs_fn, X=abs_nodes, header=str(n_abs), comments='',
           fmt=['%-12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i'])

abs_fn = os.path.join(home_dir, 'absorbing_surface_file_ymax' + name_suffix_8)
with open(abs_fn, 'r') as f_abs:
    n_abs = f_abs.readline()
n_abs = int(n_abs)
n_abs_total += n_abs
abs_nodes = np.loadtxt(fname=abs_fn, dtype=int, skiprows=1)
element_abs = abs_nodes[:, 0]
abs_nodes = np.concatenate((np.reshape(element_abs, newshape=(n_abs, 1)), mesh27_nodes[np.ix_(element_abs - 1, idx_ymax)]), axis=1)
abs_fn = os.path.join(home_dir, 'absorbing_surface_file_ymax' + name_suffix_27)
np.savetxt(fname=abs_fn, X=abs_nodes, header=str(n_abs), comments='',
           fmt=['%-12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i'])

abs_fn = os.path.join(home_dir, 'absorbing_surface_file_ymin' + name_suffix_8)
with open(abs_fn, 'r') as f_abs:
    n_abs = f_abs.readline()
n_abs = int(n_abs)
n_abs_total += n_abs
abs_nodes = np.loadtxt(fname=abs_fn, dtype=int, skiprows=1)
element_abs = abs_nodes[:, 0]
abs_nodes = np.concatenate((np.reshape(element_abs, newshape=(n_abs, 1)), mesh27_nodes[np.ix_(element_abs - 1, idx_ymin)]), axis=1)
abs_fn = os.path.join(home_dir, 'absorbing_surface_file_ymin' + name_suffix_27)
np.savetxt(fname=abs_fn, X=abs_nodes, header=str(n_abs), comments='',
           fmt=['%-12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i'])

abs_fn = os.path.join(home_dir, 'free_or_absorbing_surface_file_zmax' + name_suffix_8)
with open(abs_fn, 'r') as f_abs:
    n_abs = f_abs.readline()
n_abs = int(n_abs)
n_abs_total += n_abs
abs_nodes = np.loadtxt(fname=abs_fn, dtype=int, skiprows=1)
element_abs = abs_nodes[:, 0]
abs_nodes = np.concatenate((np.reshape(element_abs, newshape=(n_abs, 1)), mesh27_nodes[np.ix_(element_abs - 1, idx_top)]), axis=1)
abs_fn = os.path.join(home_dir, 'free_or_absorbing_surface_file_zmax' + name_suffix_27)
np.savetxt(fname=abs_fn, X=abs_nodes, header=str(n_abs), comments='',
           fmt=['%-12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i', '%12i'])

print('counted boundary faces: ' + str(len_dict))
print('actual boundary faces: ' + str(n_abs_total))
