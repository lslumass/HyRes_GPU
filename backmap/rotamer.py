import numpy as np
import math

def norm_vector(v):
    """
    Normalise a vector.
    
    Parameters
    ----------
    v : ndarray
        The array containg the vector(s).
        The vectors are represented by the last axis.
    """
    factor = np.linalg.norm(v, axis=-1)
    if isinstance(factor, np.ndarray):
        v /= factor[..., np.newaxis]
    else:
        v /= factor

def cal_normal(v1, v2):
    return np.cross(v1, v2)

def cal_distance(coord1, coord2):
    temp = np.array(coord1) - np.array(coord2)
    euclid_dist = np.sqrt(np.dot(temp.T, temp))
    return euclid_dist


def cal_angle(v1, v2):
    dotproc = sum((a*b) for a, b in zip(v1, v2))
    len1 = math.sqrt(sum((a*b) for a, b in zip(v1, v1)))
    len2 = math.sqrt(sum((a*b) for a, b in zip(v2, v2)))
    return math.acos(dotproc/(len1*len2))/np.pi*180.0
    

def rotate_about_axis(coords, axis, angle, support=None):
    """
    Rotate the given atoms or coordinates about a given axis by a given
    angle.
    
    Parameters
    ----------
    atoms : Atom or AtomArray or AtomArrayStack or ndarray, shape=(3,) or shape=(n,3) or shape=(m,n,3)
        The atoms of which the coordinates are transformed.
        The coordinates can be directly provided as :class:`ndarray`.
    axis : array-like, length=3
        A vector representing the direction of the rotation axis.
        The length of the vector is irrelevant.
    angle : float
        The rotation angle in radians.
    support : array-like, length=3, optional
        An optional support vector for the rotation axis, i.e. the
        center of the rotation.
        By default, the center of the rotation is at *(0,0,0)*.
    
    Returns
    -------
    transformed : Atom or AtomArray or AtomArrayStack or ndarray, shape=(3,) or shape=(n,3) or shape=(m,n,3)
        A copy of the input atoms or coordinates, rotated about the
        given axis.
        
    See Also
    --------
    rotate
    rotate_centered

    Examples
    --------
    Rotation about a custom axis on the *y*-*z*-plane by 90 degrees:

    >>> position = np.array([2.0, 0.0, 0.0])
    >>> axis = [0.0, 1.0, 1.0]
    >>> rotated = rotate_about_axis(position, axis, angle=0.5*np.pi)
    >>> print(rotated)
    [ 1.225e-16  1.414e+00 -1.414e+00]
    """
    positions = coords
    if support is not None:
        # Transform coordinates
        # so that the axis support vector is at (0,0,0)
        positions -= np.asarray(support)
    
    # Normalize axis
    axis = np.asarray(axis, dtype=np.float32).copy()
    if np.linalg.norm(axis) == 0:
        raise ValueError("Length of the rotation axis is 0")
    norm_vector(axis)
    # Save some interim values, that are used repeatedly in the
    # rotation matrix, for clarity and to save computation time
    sin_a = np.sin(angle)
    cos_a = np.cos(angle)
    icos_a = 1 - cos_a
    x = axis[...,0]
    y = axis[...,1]
    z = axis[...,2]
    # Rotation matrix is taken from
    # https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
    rot_matrix = np.array([
        [ cos_a + icos_a*x**2,  icos_a*x*y - z*sin_a,  icos_a*x*z + y*sin_a],
        [icos_a*x*y + z*sin_a,   cos_a + icos_a*y**2,  icos_a*y*z - x*sin_a],
        [icos_a*x*z - y*sin_a,  icos_a*y*z + x*sin_a,   cos_a + icos_a*z**2]
    ])

    # For proper rotation reshape into a maximum of 2 dimensions
    orig_ndim = positions.ndim
    if orig_ndim > 2:
        orig_shape = positions.shape
        positions = positions.reshape(-1, 3)
    # Apply rotation
    positions = np.dot(rot_matrix, positions.T).T
    # Reshape back into original shape
    if orig_ndim > 2:
        positions = positions.reshape(*orig_shape)

    if support is not None:
        # Transform coordinates back to original support vector position
        positions += np.asarray(support)
    
    return positions


def opt_side_chain(resname, refs, sides):
    """
    ref: the coords of hyres bead as the reference
    coords: the coords of atomic side chains to rotate
    """

    # for SER, THR, CYS, VAL
    if resname in ['SER', 'THR', 'CYS', 'VAL']:
        CA_r = refs.atoms[0].position
        CB_r = refs.atoms[1].position
        CA = sides.select_atoms("name CA").atoms[0].position
        CB = sides.select_atoms("name CB").atoms[0].position
        axis = np.array(CB) - np.array(CA)
        rotations = sides.select_atoms("not name N CA C O HA")

        angle = 180
        opt_positions = rotations.atoms.positions
        v1 = np.array(CB_r) - np.array(CA)
        for theta in range(0, 360, 5):
            theta = theta/180.0*np.pi
            for atom in rotations.atoms:
                atom.position = rotate_about_axis(atom.position, axis, theta, CA)

            com = rotations.center_of_mass()
            v2 = np.array(com) - np.array(CA)
            if cal_angle(v1, v2) <= angle:
                angle = cal_angle(v1, v2)
                opt_positions = rotations.atoms.positions

        rotations.atoms.positions = opt_positions

    # for ASP, ASN, LEU
    if resname in ['ASP', 'ASN', 'LEU', 'GLU', 'GLN', 'MET']:
        CA_r = refs.atoms[0].position
        CB_r = refs.atoms[1].position
        CA = sides.select_atoms("name CA").atoms[0].position
        CB = sides.select_atoms("name CB").atoms[0].position
        CA_CB = np.array(CB) - np.array(CA)
        rotations = sides.select_atoms("not name N CA C O HA")

        angle = 180
        opt_positions = rotations.atoms.positions
        v1 = np.array(CB_r) - np.array(CA)
        for theta in range(0, 360, 5):
            theta = theta/180.0*np.pi
            for atom in rotations.atoms:
                atom.position = rotate_about_axis(atom.position, CA_CB, theta, CA)
            
            CC = sides.select_atoms("name CC").atoms[0].position
            CB_CC = np.array(CC) - np.array(CB)
            part2 = sides.select_atoms("not name N CA C O HA CB HB")
            for phi in range(0, 360, 5):
                phi = phi/180.0*np.pi
                for atom in part2.atoms:
                    atom.position = rotate_about_axis(atom.position, CB_CC, phi, CB)
                
                com = rotations.center_of_mass()
                v2 = np.array(com) - np.array(CA)
                if cal_angle(v1, v2) <= angle:
                    angle = cal_angle(v1, v2)
                    steps = sides.select_atoms("not name N CA C O HA")
                    opt_positions = steps.atoms.positions

        rotations.atoms.positions = opt_positions

    # for ILE
    if resname in ['ILE']:
        CA_r = refs.atoms[0].position
        CB_r = refs.atoms[1].position
        CA = sides.select_atoms("name CA").atoms[0].position
        CB = sides.select_atoms("name CB").atoms[0].position
        CA_CB = np.array(CB) - np.array(CA)
        rotations = sides.select_atoms("not name N CA C O HA")

        angle = 180
        opt_positions = rotations.atoms.positions
        v1 = np.array(CB_r) - np.array(CA)
        for theta in range(0, 360, 5):
            theta = theta/180.0*np.pi
            for atom in rotations.atoms:
                atom.position = rotate_about_axis(atom.position, CA_CB, theta, CA)
            
            CC = sides.select_atoms("name CC").atoms[0].position
            CB_CC = np.array(CC) - np.array(CB)
            part2 = sides.select_atoms("name CF HF")
            for phi in range(0, 360, 5):
                phi = phi/180.0*np.pi
                for atom in part2.atoms:
                    atom.position = rotate_about_axis(atom.position, CB_CC, phi, CB)
                
                com = rotations.center_of_mass()
                v2 = np.array(com) - np.array(CA)
                if cal_angle(v1, v2) <= angle:
                    angle = cal_angle(v1, v2)
                    steps = sides.select_atoms("not name N CA C O HA")
                    opt_positions = steps.atoms.positions

        rotations.atoms.positions = opt_positions

    # for ARG
    if resname in ['ARG']:
        CA_r = refs.select_atoms("name CA").atoms[0].position
        CB_r = refs.select_atoms("name CB").atoms[0].position
        CC_r = refs.select_atoms("name CC").atoms[0].position
        CA = sides.select_atoms("name CA").atoms[0].position
        CB = sides.select_atoms("name CB").atoms[0].position
        CA_CB = np.array(CB) - np.array(CA)
        rotations = sides.select_atoms("not name N CA C O HA")

        angle = 180
        opt_positions = rotations.atoms.positions
        v1 = np.array(CB_r) - np.array(CA)
        for theta in range(0, 360, 5):
            theta = theta/180.0*np.pi
            for atom in rotations.atoms:
                atom.position = rotate_about_axis(atom.position, CA_CB, theta, CA)
            
            CC = sides.select_atoms("name CC").atoms[0].position
            CB_CC = np.array(CC) - np.array(CB)
            part2 = sides.select_atoms("not name N CA C O HA CB HB")
            for phi in range(0, 360, 5):
                phi = phi/180.0*np.pi
                for atom in part2.atoms:
                    atom.position = rotate_about_axis(atom.position, CB_CC, phi, CB)
                
                grp_CB = sides.select_atoms("name CB HB CC HC CD HD")
                com = rotations.center_of_mass()
                v2 = np.array(com) - np.array(CA)
                if cal_angle(v1, v2) <= angle:
                    angle = cal_angle(v1, v2)
                    steps = sides.select_atoms("not name N CA C O HA")
                    opt_positions = steps.atoms.positions
        rotations.atoms.positions = opt_positions

        v1 = np.array(CC_r) - np.array(CA)
        CD = sides.select_atoms("name CD").atoms[0].position
        N1 = sides.select_atoms("name N1").atoms[0].position
        CD_N1 = np.array(N1) - np.array(CD)
        part3 = sides.select_atoms("not name N CA C O HA CB HB CC HC CD HD")
        opt_positions = rotations.atoms.positions
        angle = 180
        for theta in range(0, 360, 5):
            theta = theta/180.0*np.pi
            for atom in part3.atoms:
                atom.position = rotate_about_axis(atom.position, CD_N1, theta, CD)
            
            com = part3.center_of_mass()
            v2 = np.array(com) - np.array(CA)
            if cal_angle(v1, v2) <= angle:
                angle = cal_angle(v1, v2)
                opt_positions = part3.atoms.positions
        part3.atoms.positions = opt_positions

    # for LYS
    if resname in ['LYS']:
        CA_r = refs.select_atoms("name CA").atoms[0].position
        CB_r = refs.select_atoms("name CB").atoms[0].position
        CC_r = refs.select_atoms("name CC").atoms[0].position
        CA = sides.select_atoms("name CA").atoms[0].position
        CB = sides.select_atoms("name CB").atoms[0].position
        CA_CB = np.array(CB) - np.array(CA)
        rotations = sides.select_atoms("not name N CA C O HA")

        angle = 180
        opt_positions = rotations.atoms.positions
        v1 = np.array(CB_r) - np.array(CA)
        for theta in range(0, 360, 5):
            theta = theta/180.0*np.pi
            for atom in rotations.atoms:
                atom.position = rotate_about_axis(atom.position, CA_CB, theta, CA)
            
            CC = sides.select_atoms("name CC").atoms[0].position
            CB_CC = np.array(CC) - np.array(CB)
            part2 = sides.select_atoms("not name N CA C O HA CB HB")
            for phi in range(0, 360, 5):
                phi = phi/180.0*np.pi
                for atom in part2.atoms:
                    atom.position = rotate_about_axis(atom.position, CB_CC, phi, CB)
                
                grp_CB = sides.select_atoms("name CB HB CC HC CD HD")
                com = rotations.center_of_mass()
                v2 = np.array(com) - np.array(CA)
                if cal_angle(v1, v2) <= angle:
                    angle = cal_angle(v1, v2)
                    steps = sides.select_atoms("not name N CA C O HA")
                    opt_positions = steps.atoms.positions
        rotations.atoms.positions = opt_positions

        v1 = np.array(CC_r) - np.array(CA)
        CC = sides.select_atoms("name CC").atoms[0].position
        CD = sides.select_atoms("name CD").atoms[0].position
        CC_CD = np.array(CD) - np.array(CC)
        part3 = sides.select_atoms("name CD HD CE HE N1 HN")
        opt_positions = part3.atoms.positions
        angle = 180
        for theta in range(0, 360, 5):
            theta = theta/180.0*np.pi
            for atom in part3.atoms:
                atom.position = rotate_about_axis(atom.position, CC_CD, theta, CC)
            
            grp_CC = sides.select_atoms("name CE HE N1 HN")
            com = grp_CC.center_of_mass()
            v2 = np.array(com) - np.array(CA)
            if cal_angle(v1, v2) <= angle:
                angle = cal_angle(v1, v2)
                opt_positions = part3.atoms.positions
        part3.atoms.positions = opt_positions
        

# reserved for high_precise MET
#    # for MET
#    if resname in ['MET']:
#        CA_r = refs.atoms[0].position
#        CB_r = refs.atoms[1].position
#        CA = sides.select_atoms("name CA").atoms[0].position
#        CB = sides.select_atoms("name CB").atoms[0].position
#        CA_CB = np.array(CB) - np.array(CA)
#        rotations = sides.select_atoms("not name N CA C O HA")
#
#        angle = 180
#        opt_positions = rotations.atoms.positions
#        v1 = np.array(CB_r) - np.array(CA)
#        for theta in range(0, 360, 5):
#            theta = theta/180.0*np.pi
#            for atom in rotations.atoms:
#                atom.position = rotate_about_axis(atom.position, CA_CB, theta, CA)
#            
#            CC = sides.select_atoms("name CC").atoms[0].position
#            CB_CC = np.array(CC) - np.array(CB)
#            part2 = sides.select_atoms("not name N CA C O HA CB HB")
#            for phi in range(0, 360, 5):
#                phi = phi/180.0*np.pi
#                for atom in part2.atoms:
#                    atom.position = rotate_about_axis(atom.position, CB_CC, phi, CB)
#                
#                S = sides.select_atoms("name S").atoms[0].position
#                CC_S = np.array(S) - np.array(CC)
#                part3 = sides.select_atoms("not name N CA C O HA CB HB CC HC")
#                for psi in range(0, 360, 5):
#                    psi = psi/180.0*np.pi
#                    for atom in part3.atoms:
#                        atom.position = rotate_about_axis(atom.position, CC_S, psi, CC)
#                
#                    com = rotations.center_of_mass()
#                    v2 = np.array(com) - np.array(CA)
#                    if cal_angle(v1, v2) <= angle:
#                        angle = cal_angle(v1, v2)
#                        steps = sides.select_atoms("not name N CA C O HA")
#                        opt_positions = steps.atoms.positions
#
#        rotations.atoms.positions = opt_positions

    # for HIS, PHE, TYR
    if resname in ['HIS', 'PHE', 'TYR', 'TRP']:
        CA_r = refs.select_atoms("name CA").atoms[0].position
        CB_r = refs.select_atoms("name CB").atoms[0].position
        CC_r = refs.select_atoms("name CC").atoms[0].position
        CD_r = refs.select_atoms("name CD").atoms[0].position
        CA = sides.select_atoms("name CA").atoms[0].position
        CB = sides.select_atoms("name CB").atoms[0].position
        CA_CB = np.array(CB) - np.array(CA)
        rotations = sides.select_atoms("not name N CA C O HA")

        angle = 180
        opt_positions = rotations.atoms.positions
        v1 = np.array(CB_r) - np.array(CA)
        for theta in range(0, 360, 5):
            theta = theta/180.0*np.pi
            for atom in rotations.atoms:
                atom.position = rotate_about_axis(atom.position, CA_CB, theta, CA)
            
            grp_CB = sides.select_atoms("name CB HB CC")
            com_CB = grp_CB.center_of_mass()
            v2 = np.array(com_CB) - np.array(CA)
            if cal_angle(v1, v2) <= angle:
                angle = cal_angle(v1, v2)
                opt_positions = rotations.atoms.positions
        rotations.atoms.positions = opt_positions

        CC = sides.select_atoms("name CC").atoms[0].position
        CB_CC = np.array(CC) - np.array(CB)
        part2 = sides.select_atoms("not name N CA C O HA CB HB")
        cut_off = 180.0
        nomal1 = cal_normal(np.array(CC_r)-np.array(CB_r), np.array(CD_r)-np.array(CB_r))
        for phi in range(0, 360, 5):
            phi = phi/180.0*np.pi
            for atom in part2.atoms:
                atom.position = rotate_about_axis(atom.position, CB_CC, phi, CB)
            
            if resname == 'HIS':
                grp_CB = sides.select_atoms("name CB HB CC")
                grp_CC = sides.select_atoms("name N1 CE HE")
                grp_CD = sides.select_atoms("name CD HD N2 HN")
            elif resname == 'PHE':
                grp_CB = sides.select_atoms("name CB HB CC")
                grp_CC = sides.select_atoms("name CE HE CG HG")
                grp_CD = sides.select_atoms("name CF HF CH HH")
            elif resname == 'TYR':
                grp_CB = sides.select_atoms("name CB HB CC")
                grp_CC = sides.select_atoms("name CE HE CG HG")
                grp_CD = sides.select_atoms("name CF HF CH O1 HO")
            elif resname == 'TRP':
                grp_CB = sides.select_atoms("name CB HB CC")
                grp_CC = sides.select_atoms("name CD HD N1 HN")
                grp_CD = sides.select_atoms("name CE CF")
                grp_CE = sides.select_atoms("name CH HH CJ HJ")
                grp_CF = sides.select_atoms("name CG HG CI HI")
            
            com_B = grp_CB.center_of_mass()
            com_C = grp_CC.center_of_mass()
            com_D = grp_CD.center_of_mass()
            distance = cal_distance(com_B, CB_r) + cal_distance(com_C, CC_r) + cal_distance(com_D, CD_r)
            if distance <= cut_off:
                cut_off = distance
                opt_positions = part2.atoms.positions

        part2.atoms.positions = opt_positions

