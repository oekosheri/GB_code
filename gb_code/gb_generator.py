
# !/usr/bin/env python
"""
This module produces GB structures. You need to run csl_generator first
to get the info necessary for your grain boundary of interest.

 The code runs in one mode only and takes all the necessary options to write
 a final GB structure from the input io_file, which is written  after you run
 csl_generator.py. You must customize the io_file that comes with some default
 values. For ex.: input the GB_plane of interest from  running the
 CSl_generator in the second mode. Once you have completed customizing the
 io_file, run:

 'gb_generator.py io_file'
 """

import sys
import numpy as np
from numpy import dot, cross
from numpy.linalg import det, norm
import csl_generator as cslgen
import warnings


class GB_character:
    """
    a grain boundary class, encompassing all the characteristics of a GB.
    """
    def __init__(self):
        self.axis = np.array([1, 0, 0])
        self.sigma = 1
        self.theta = 0
        self.m = 1
        self.n = 1
        self.R = np.eye(1)
        self.basis = 'fcc'
        self.LatP = 4.05
        self.gbplane = np.array([1, 1, 1])
        self.ortho1 = np.eye(3)
        self.ortho2 = np.eye(3)
        self.ortho = np.eye(3)
        self.atoms = np.eye(3)
        self.atoms1 = np.eye(3)
        self.atoms2 = np.eye(3)
        self.rot1 = np.eye(3)
        self.rot2 = np.eye(3)
        self.Num = 0
        self.dim = np.array([1, 1, 1])
        self.overD = 0
        self.whichG = 'g1'
        self.trans = False
        self.File = 'LAMMPS'

    def ParseGB(self, axis, basis, LatP, m, n, gb):
        """
        parses the GB input axis, basis, lattice parameter,
        m and n integers and gb plane.
        """
        self.axis = np.array(axis)
        self.m = int(m)
        self.n = int(n)
        self.sigma = cslgen.get_cubic_sigma(self.axis, self.m, self.n)
        self.theta = cslgen.get_cubic_theta(self.axis, self.m, self.n)
        self.R = cslgen.rot(self.axis, self.theta)

        if (str(basis) == 'fcc' or str(basis) == 'bcc' or str(basis) == 'sc' or
           str(basis) == 'diamond'):

            self.basis = str(basis)

            self.LatP = float(LatP)
            self.gbplane = np.array(gb)

            try:
                self.ortho1, self.ortho2, self.Num = \
                    cslgen.Find_Orthogonal_cell(self.basis,
                                                self.axis,
                                                self.m,
                                                self.n,
                                                self.gbplane)

            except:
                print("""
                    Could not find the orthogonal cells.... Most likely the
                    input GB_plane is "NOT" a CSL plane. Go back to the first
                    script and double check!
                    """)
                sys.exit()
        else:
            print("Sorry! For now only works for cubic lattices ... ")
            sys.exit()

    def WriteGB(self, overlap=0.0, rigid=False,
                dim1=1, dim2=1, dim3=1, file='LAMMPS',
                **kwargs):
        """
        parses the arguments overlap distance, dimensions, rigid body translate
        and file type and writes the final structure to the file.
        Possible keys:
            (whichG, a, b)
        """
        self.overD = float(overlap)
        self.trans = rigid
        self.dim = np.array([int(dim1), int(dim2), int(dim3)])
        self.File = file
        if self.overD > 0:
            try:
                self.whichG = kwargs['whichG']
            except:
                print('decide on whichG!')
                sys.exit()
            if self.trans:
                try:
                    a = int(kwargs['a'])
                    b = int(kwargs['b'])
                except:
                    print('Make sure the a and b integers are there!')
                    sys.exit()
            self.Expand_Super_cell()
            xdel, _, x_indice, y_indice = self.Find_overlapping_Atoms()
            print("<<------ {} atoms are being removed! ------>>"
                  .format(len(xdel)))

            if self.whichG == "G1" or self.whichG == "g1":
                self.atoms1 = np.delete(self.atoms1, x_indice, axis=0)

            elif self.whichG == "G2" or self.whichG == "g2":
                self.atoms2 = np.delete(self.atoms2, y_indice, axis=0)

            else:
                print("You must choose either 'g1', 'g2' ")
                sys.exit()
            if not self.trans:
                count = 0
                print("<<------ 1 GB structure is being created! ------>>")
                if self.File == "LAMMPS":
                    self.Write_to_Lammps(count)
                elif self.File == "VASP":
                    self.Write_to_Vasp(count)
                else:
                    print("The output file must be either LAMMPS or VASP!")
            elif self.trans:
                self.Translate(a, b)

        elif self.overD == 0:
            if self.trans:
                try:
                    a = int(kwargs['a'])
                    b = int(kwargs['b'])
                except:
                    print('Make sure the a and b integers are there!')
                    sys.exit()
                print("<<------ 0 atoms are being removed! ------>>")
                self.Expand_Super_cell()
                self.Translate(a, b)

            else:
                self.Expand_Super_cell()
                count = 0
                print("<<------ 1 GB structure is being created! ------>>")
                if self.File == "LAMMPS":
                    self.Write_to_Lammps(count)
                elif self.File == "VASP":
                    self.Write_to_Vasp(count)
                else:
                    print("The output file must be either LAMMPS or VASP!")
        else:
            print('Overlap distance is not inputted incorrectly!')
            sys.exit()

    def CSL_Ortho_unitcell_atom_generator(self):

        """
        populates a unitcell from the orthogonal vectors.
        """
        Or = self.ortho.T
        Orint = cslgen.integerMatrix(Or)
        LoopBound = np.zeros((3, 2), dtype=float)
        transformed = []
        CubeCoords = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0],
                              [0, 1, 1], [1, 0, 1], [1, 1, 1], [0, 0, 0]],
                              dtype=float)
        for i in range(len(CubeCoords)):
            transformed.append(np.dot(Orint.T, CubeCoords[i]))

        # Finding bounds for atoms in a CSL unitcell:

        LoopBound[0, :] = [min(np.array(transformed)[:, 0]),
                           max(np.array(transformed)[:, 0])]
        LoopBound[1, :] = [min(np.array(transformed)[:, 1]),
                           max(np.array(transformed)[:, 1])]
        LoopBound[2, :] = [min(np.array(transformed)[:, 2]),
                           max(np.array(transformed)[:, 2])]

        # Filling up the unitcell:

        Tol = 1
        x = np.arange(LoopBound[0, 0] - Tol, LoopBound[0, 1] + Tol + 1, 1)
        y = np.arange(LoopBound[1, 0] - Tol, LoopBound[1, 1] + Tol + 1, 1)
        z = np.arange(LoopBound[2, 0] - Tol, LoopBound[2, 1] + Tol + 1, 1)
        V = len(x) * len(y) * len(z)
        indice = (np.stack(np.meshgrid(x, y, z)).T).reshape(V, 3)
        Base = cslgen.Basis(str(self.basis))
        Atoms = []
        tol = 0.001
        if V > 5e6:
            print("Warning! It may take a very long time"
                  "to produce this cell!")
        # produce Atoms:

        for i in range(V):
            for j in range(len(Base)):
                Atoms.append(indice[i, 0:3] + Base[j, 0:3])
        Atoms = np.array(Atoms)

        # Cell conditions
        Con1 = dot(Atoms, Or[0]) / norm(Or[0]) + tol
        Con2 = dot(Atoms, Or[1]) / norm(Or[1]) + tol
        Con3 = dot(Atoms, Or[2]) / norm(Or[2]) + tol
        # Application of the conditions:
        Atoms = (Atoms[(Con1 >= 0) & (Con1 <= norm(Or[0])) & (Con2 >= 0) &
                 (Con2 <= norm(Or[1])) &
                 (Con3 >= 0) & (Con3 <= norm(Or[2]))])

        if len(Atoms) == (round(det(Or) * len(Base), 7)).astype(int):
            self.Atoms = Atoms
        else:
            self.Atoms = None
        return

    def CSL_Bicrystal_Atom_generator(self):
        """
        builds the unitcells for both grains g1 and g2.
        """
        Or_1 = self.ortho1.T
        Or_2 = self.ortho2.T
        self.rot1 = np.array([Or_1[0, :] / norm(Or_1[0, :]),
                             Or_1[1, :] / norm(Or_1[1, :]),
                             Or_1[2, :] / norm(Or_1[2, :])])
        self.rot2 = np.array([Or_2[0, :] / norm(Or_2[0, :]),
                             Or_2[1, :] / norm(Or_2[1, :]),
                             Or_2[2, :] / norm(Or_2[2, :])])

        self.ortho = self.ortho1.copy()
        self.CSL_Ortho_unitcell_atom_generator()
        self.atoms1 = self.Atoms

        self.ortho = self.ortho2.copy()
        self.CSL_Ortho_unitcell_atom_generator()
        self.atoms2 = self.Atoms

        self.atoms1 = dot(self.rot1, self.atoms1.T).T
        self.atoms2 = dot(self.rot2, self.atoms2.T).T
        # tol = 0.01
        self.atoms2[:, 0] = self.atoms2[:, 0] - norm(Or_2[0, :])  # - tol
        # print(self.atoms2, norm(Or_2[0, :]) )
        return

    def Expand_Super_cell(self):
        """
        expands the smallest CSL unitcell to the given dimensions.
        """
        a = norm(self.ortho1[:, 0])
        b = norm(self.ortho1[:, 1])
        c = norm(self.ortho1[:, 2])
        dimX, dimY, dimZ = self.dim

        X = self.atoms1.copy()
        Y = self.atoms2.copy()

        X_new = []
        Y_new = []
        for i in range(dimX):
            for j in range(dimY):
                for k in range(dimZ):
                    Position1 = [i * a, j * b, k * c]
                    Position2 = [-i * a, j * b, k * c]
                    for l in range(len(X)):
                        X_new.append(Position1[0:3] + X[l, 0:3])
                    for m in range(len(Y)):
                        Y_new.append(Position2[0:3] + Y[m, 0:3])

        self.atoms1 = np.array(X_new)
        self.atoms2 = np.array(Y_new)

        return

    def Find_overlapping_Atoms(self):

        """
        finds the overlapping atoms.
        """
        periodic_length = norm(self.ortho1[:, 0]) * self.dim[0]
        periodic_image = self.atoms2 + [periodic_length * 2, 0, 0]
        # select atoms contained in a smaller window around the GB and its
        # periodic image
        IndX = np.where([(self.atoms1[:, 0] < 1) |
                         (self.atoms1[:, 0] > (periodic_length - 1))])[1]
        IndY = np.where([self.atoms2[:, 0] > -1])[1]
        IndY_image = np.where([periodic_image[:, 0] <
                              (periodic_length + 1)])[1]
        X_new = self.atoms1[IndX]
        Y_new = np.concatenate((self.atoms2[IndY], periodic_image[IndY_image]))
        IndY_new = np.concatenate((IndY, IndY_image))
        # create a meshgrid search routine
        x = np.arange(0, len(X_new), 1)
        y = np.arange(0, len(Y_new), 1)
        indice = (np.stack(np.meshgrid(x, y)).T).reshape(len(x) * len(y), 2)
        norms = norm(X_new[indice[:, 0]] - Y_new[indice[:, 1]], axis=1)
        indice_x = indice[norms < self.overD][:, 0]
        indice_y = indice[norms < self.overD][:, 1]
        X_del = X_new[indice_x]
        Y_del = Y_new[indice_y]

        if (len(X_del) != len(Y_del)):
            print("Warning! the number of deleted atoms"
                  "in the two grains are not equal!")
        # print(type(IndX), len(IndY), len(IndY_image))
        return (X_del, Y_del, IndX[indice_x], IndY_new[indice_y])

    def Translate(self, a, b):

        """
        translates the GB on a mesh created by a, b integers and writes
        to LAMMPS or VASP.
        """
        tol = 0.001
        if (1 - cslgen.ang(self.gbplane, self.axis) < tol):

            M1, _ = cslgen.Create_minimal_cell_Method_1(
                     self.sigma, self.axis, self.R)
            D = (1 / self.sigma * cslgen.DSC_vec(self.basis, self.sigma, M1))
            Dvecs = cslgen.DSC_on_plane(D, self.gbplane)
            TransDvecs = np.round(dot(self.rot1, Dvecs), 7)
            shift1 = TransDvecs[:, 0] / 2
            shift2 = TransDvecs[:, 1] / 2
            a = b = 3
        else:
            # a = 10
            # b = 5
            if norm(self.ortho1[:, 1]) > norm(self.ortho1[:, 2]):

                shift1 = (1 / a) * (norm(self.ortho1[:, 1]) *
                                    np.array([0, 1, 0]))
                shift2 = (1 / b) * (norm(self.ortho1[:, 2]) *
                                    np.array([0, 0, 1]))
            else:
                shift1 = (1 / a) * (norm(self.ortho1[:, 2]) *
                                    np.array([0, 0, 1]))
                shift2 = (1 / b) * (norm(self.ortho1[:, 1]) *
                                    np.array([0, 1, 0]))
        print("<<------ {} GB structures are being created! ------>>"
              .format(int(a*b)))

        XX = self.atoms1
        count = 0
        if self.File == 'LAMMPS':

            for i in range(a):
                for j in range(b):
                    count += 1
                    shift = i * shift1 + j * shift2
                    atoms1_new = XX.copy() + shift
                    self.atoms1 = atoms1_new
                    self.Write_to_Lammps(count)
        elif self.File == 'VASP':

            for i in range(a):
                for j in range(b):
                    count += 1
                    shift = i * shift1 + j * shift2
                    atoms1_new = XX.copy() + shift
                    self.atoms1 = atoms1_new
                    self.Write_to_Vasp(count)
        else:
            print("The output file must be either LAMMPS or VASP!")

    def Write_to_Vasp(self, trans):
        """
        write a single GB without translations to POSCAR.
        """
        name = 'POS_G'
        plane = str(self.gbplane[0])+str(self.gbplane[1])+str(self.gbplane[2])
        if self.overD > 0:
            overD = str(self.overD)
        else:
            overD = str(None)
        Trans = str(trans)
        # tol = 0.001
        X = self.atoms1.copy()
        Y = self.atoms2.copy()
        X_new = X * self.LatP
        Y_new = Y * self.LatP
        dimx, dimy, dimz = self.dim

        xlo = -1 * np.round(norm(self.ortho1[:, 0]) * dimx * self.LatP, 8)
        xhi = np.round(norm(self.ortho1[:, 0]) * dimx * self.LatP, 8)
        LenX = xhi - xlo
        ylo = 0.0
        yhi = np.round(norm(self.ortho1[:, 1]) * dimy * self.LatP, 8)
        LenY = yhi - ylo
        zlo = 0.0
        zhi = np.round(norm(self.ortho1[:, 2]) * dimz * self.LatP, 8)
        LenZ = zhi - zlo

        Wf = np.concatenate((X_new, Y_new))

        with open(name + plane + '_' + overD + '_' + Trans, 'w') as f:
            f.write('#POSCAR written by GB_code \n')
            f.write('1 \n')
            f.write('{0:.8f} 0.0 0.0 \n'.format(LenX))
            f.write('0.0 {0:.8f} 0.0 \n'.format(LenY))
            f.write('0.0 0.0 {0:.8f} \n'.format(LenZ))
            f.write('{} {} \n'.format(len(X), len(Y)))
            f.write('Cartesian\n')
            np.savetxt(f, Wf, fmt='%.8f %.8f %.8f')
        f.close()

    def Write_to_Lammps(self, trans):
        """
        write a single GB without translations to LAMMPS.
        """
        name = 'input_G'
        plane = str(self.gbplane[0])+str(self.gbplane[1])+str(self.gbplane[2])
        if self.overD > 0:
            overD = str(self.overD)
        else:
            overD = str(None)
        Trans = str(trans)
        # tol = 0.001
        X = self.atoms1.copy()
        Y = self.atoms2.copy()

        NumberAt = len(X) + len(Y)
        X_new = X * self.LatP
        Y_new = Y * self.LatP
        dimx, dimy, dimz = self.dim

        xlo = -1 * np.round(norm(self.ortho1[:, 0]) * dimx * self.LatP, 8)
        xhi = np.round(norm(self.ortho1[:, 0]) * dimx * self.LatP, 8)
        ylo = 0.0
        yhi = np.round(norm(self.ortho1[:, 1]) * dimy * self.LatP, 8)
        zlo = 0.0
        zhi = np.round(norm(self.ortho1[:, 2]) * dimz * self.LatP, 8)

        Type1 = np.ones(len(X_new), int).reshape(1, -1)
        Type2 = 2 * np.ones(len(Y_new), int).reshape(1, -1)
        # Type = np.concatenate((Type1, Type2), axis=1)

        Counter = np.arange(1, NumberAt + 1).reshape(1, -1)

        # data = np.concatenate((X_new, Y_new))
        W1 = np.concatenate((Type1.T, X_new), axis=1)
        W2 = np.concatenate((Type2.T, Y_new), axis=1)
        Wf = np.concatenate((W1, W2))
        FinalMat = np.concatenate((Counter.T, Wf), axis=1)

        with open(name + plane + '_' + overD + '_' + Trans, 'w') as f:
            f.write('#Header \n \n')
            f.write('{} atoms \n \n'.format(NumberAt))
            f.write('2 atom types \n \n')
            f.write('{0:.8f} {1:.8f} xlo xhi \n'.format(xlo, xhi))
            f.write('{0:.8f} {1:.8f} ylo yhi \n'.format(ylo, yhi))
            f.write('{0:.8f} {1:.8f} zlo zhi \n\n'.format(zlo, zhi))
            f.write('Atoms \n \n')
            np.savetxt(f, FinalMat, fmt='%i %i %.8f %.8f %.8f')
        f.close()

    def __str__(self):
        return "GB_character"


def main():
    import yaml
    if len(sys.argv) == 2:
        io_file = sys.argv[1]
        file = open(io_file, 'r')
        in_params = yaml.load(file)

        try:
            axis = np.array(in_params['axis'])
            m = int(in_params['m'])
            n = int(in_params['n'])
            basis = str(in_params['basis'])
            gbplane = np.array(in_params['GB_plane'])
            LatP = in_params['lattice_parameter']
            overlap = in_params['overlap_distance']
            whichG = in_params['which_g']
            rigid = in_params['rigid_trans']
            a = in_params['a']
            b = in_params['b']
            dim1, dim2, dim3 = in_params['dimensions']
            file = in_params['File_type']

        except:
            print('Make sure the input argumnets in io_file are'
                  'put in correctly!')
            sys.exit()

        ###################

        gbI = GB_character()
        gbI.ParseGB(axis, basis, LatP, m, n, gbplane)
        gbI.CSL_Bicrystal_Atom_generator()

        if overlap > 0 and rigid:
            gbI.WriteGB(
                overlap=overlap, whichG=whichG, rigid=rigid, a=a,
                b=b, dim1=dim1, dim2=dim2, dim3=dim3, file=file
                )
        elif overlap > 0 and not rigid:
            gbI.WriteGB(
                overlap=overlap, whichG=whichG, rigid=rigid,
                dim1=dim1, dim2=dim2, dim3=dim3, file=file
                )
        elif overlap == 0 and rigid:
            gbI.WriteGB(
                overlap=overlap, rigid=rigid, a=a,
                b=b, dim1=dim1, dim2=dim2, dim3=dim3,
                file=file
                )
        elif overlap == 0 and not rigid:
            gbI.WriteGB(
                overlap=overlap, rigid=rigid,
                dim1=dim1, dim2=dim2, dim3=dim3, file=file
                )
    else:
        print(__doc__)
    return


if __name__ == '__main__':
    main()
