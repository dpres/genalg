import numpy as np
import random as random


class Atom:
    """Atom class contains atom-type objects that form clusters."""
    index = 0

    def __init__(self, pos_vector, identity):
        """Creates an atom with position coordinates and a type"""
        self._pos_vect = pos_vector
        self._identity = identity
        Atom.index += 1

    def report_pos(self):
        """Gives the coordinates of the atom"""
        return self._pos_vect

    def report_type(self):
        """Gives the type of the atom"""
        return self._identity


class Cluster:
    """Cluster class contains cluster-type objects comprising atom-type objects."""
    index = 0

    def __init__(self, list_of_atoms, dict_of_atoms):
        """Creates cluster from a list of atoms"""
        self._list_of_atoms = list_of_atoms
        self._dict_of_atoms = dict_of_atoms
        Cluster.index += 1

    def atoms_list_in_cluster(self):
        """Returns list of atoms in the cluster"""
        return self._list_of_atoms

    def atoms_dict_in_cluster(self):
        """Returns dictionary of atoms in the cluster"""
        return self._dict_of_atoms

    def cluster_energy(self):
        """Calculates the energy of a cluster by looping through all pairs of atoms in the cluster"""
        l = Cluster.atoms_list_in_cluster(self)
        d = Cluster.atoms_dict_in_cluster(self)

        energy = energy_calculator(l, d)
        return energy


def energy_calculator(atomlist, atomdict):
    """Calculates energy of a set of atoms"""
    l = atomlist
    d = atomdict

    sigma_dict = {(1, 1): sigmaAA, (1, 2): sigmaAB, (2, 1): sigmaAB, (2, 2): sigmaBB}
    energy = 0  # initialises energy
    for i in xrange(N - 1):
        for j in xrange(i + 1, N):
            type_i = Atom.report_type(d[l[i]])
            type_j = Atom.report_type(d[l[j]])
            atom_pair = (type_i, type_j)
            sigma = sigma_dict[atom_pair]  # picks appropriate sigma for the calculation

            pos_i = Atom.report_pos(d[l[i]])
            pos_j = Atom.report_pos(d[l[j]])
            d_pos = pos_j - pos_i
            r_ij = np.sqrt(np.dot(d_pos, d_pos))  # distance between centres of atoms

            frac6 = (sigma / r_ij) ** 6  # makes calculation faster

            energy += 4 * ((frac6 ** 2) - frac6)
    return energy

def mating(p1_alist, p1_adict, p2_alist, p2_adict):
    """mating step with position mutations"""
    list_of_atoms = []
    dict_of_atoms = {}

    for x in xrange(N):
        c = random.random()
        if c < 0.5:
            pos_arr = Atom.report_pos(p1_adict[p1_alist[x]]) + 0.1 * (random.random() - 0.5)
        else:
            pos_arr = Atom.report_pos(p2_adict[p2_alist[x]]) + 0.1 * (random.random() - 0.5)

        d = random.random()
        if d < 0.5:
            atom_type = Atom.report_type(p1_adict[p1_alist[x]])
        else:
            atom_type = Atom.report_type(p2_adict[p2_alist[x]])

        dict_of_atoms['Atom' + str(Atom.index)] = Atom(pos_arr, atom_type)
        list_of_atoms.append('Atom' + str(Atom.index))

    return (list_of_atoms, dict_of_atoms)








global N  # N is the total number of atoms in the cluster
global sigmaAA, sigmaAB, sigmaBB  # CoM-CoM distance of atoms XX

N = 2  # total number of atoms in the cluster
C = 4000  # starting number of clusters

# lists sigmas; sigmaBB needs fixing to allow list of sigmaBBs
sigmaAA = 1
sigmaBB = 1
sigmaAB = (sigmaAA + sigmaBB) / 2

# the size of the 2D space:
x_size = 20
y_size = 20

dict_of_clusters = {}
list_of_clusters = []  # initialises list of clusters that will get updated through the GA

while len(list_of_clusters) < C:  # creating C clusters
    atom_dict = {}
    atom_list = []

    while len(atom_list) < N:  # each cluster has N atoms
        # randomly assigns position coordinates:
        pos_arr = np.array([x_size * (random.random()), y_size * (random.random())])
        # equal probability of A/B type assignment
        # A==1, B==2
        t = random.random()
        if t < 0.5:
            atom_type = 1
        else:
            atom_type = 2
        # Atom-type object creation:
        atom_dict['Atom' + str(Atom.index)] = Atom(pos_arr, atom_type)
        atom_list.append('Atom' + str(Atom.index))
    # print atom_list
    # print Atom.report_pos(atom_dict[atom_list[0]])
    # print Atom.report_pos(atom_dict[atom_list[1]])
    # print Atom.report_type(atom_dict[atom_list[1]])

    # checks if the cluster is garbage; if not, creates cluster
    energy = energy_calculator(atom_list, atom_dict)
    if energy > 0.05:
        for atom in atom_list:
            del atom_dict[atom]  # cleans up unneeded memory
        Atom.index -= N  # rewinds the Atom index - initialisation of large clusters could cause integer overrun
    else:
        # print atom_list
        dict_of_clusters['Cluster' + str(Cluster.index)] = Cluster(atom_list, atom_dict)
        list_of_clusters.append('Cluster' + str(Cluster.index))
# print list_of_clusters
# print dict_of_clusters

#  for c in list_of_clusters:
#    print Cluster.cluster_energy(dict_of_clusters[c])

#  make a dictionary of energies of all clusters
dict_of_cluster_energies = {}
for cluster in list_of_clusters:
    dict_of_cluster_energies[cluster] = Cluster.cluster_energy(dict_of_clusters[cluster])


fail_count = 0
steps = 0
while fail_count < 10000:
    # choice of two random integers for list of clusters
    a = int(C * random.random())
    b = int(C * random.random())

    if a == b:  # try again if the same cluster is picked twice - parthenogenesis is useless
        continue

    # calls parent clusters and retrieves their atom lists/dicts ... better than retrieving them always
    parent1 = list_of_clusters[a]
    p1_alist = Cluster.atoms_list_in_cluster(dict_of_clusters[parent1])
    p1_adict = Cluster.atoms_dict_in_cluster(dict_of_clusters[parent1])

    parent2 = list_of_clusters[b]
    p2_alist = Cluster.atoms_list_in_cluster(dict_of_clusters[parent2])
    p2_adict = Cluster.atoms_dict_in_cluster(dict_of_clusters[parent2])

    # two parents produce two kids:
    for y in range(2):
        steps += 1
        (thelist, thedict) = mating(p1_alist, p1_adict, p2_alist, p2_adict)

        energy = energy_calculator(thelist, thedict)
        if energy < max(dict_of_cluster_energies.values()):
            thekey = max(dict_of_cluster_energies, key=dict_of_cluster_energies.get)
            del dict_of_clusters[thekey]
            del dict_of_cluster_energies[thekey]
            list_of_clusters.remove(thekey)
            dict_of_clusters['Cluster' + str(Cluster.index)] = Cluster(thelist, thedict)
            list_of_clusters.append('Cluster' + str(Cluster.index))
            dict_of_cluster_energies['Cluster' + str(Cluster.index)] = Cluster.cluster_energy(dict_of_clusters['Cluster' + str(Cluster.index)])
            fail_count = 0
        else:
            for atom in thelist:
                del thedict[atom]  # cleans up unneeded memory
            Atom.index -= N
            fail_count += 1
            print 'f', steps, fail_count, max(dict_of_cluster_energies.values())


print list_of_clusters
print dict_of_clusters

print dict_of_cluster_energies
