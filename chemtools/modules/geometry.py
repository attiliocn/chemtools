import numpy as np
from modules.kabsch import kabsch_algorithm

def measure_distance(p0,p1):
    b0 = p1 - p0
    distance = np.linalg.norm(b0)
    return distance

def measure_angle(p0,p1,p2):
    v1 = p0 - p1
    v2 = p2 - p1
    inner = np.inner(v1, v2)
    norms = np.linalg.norm(v1) * np.linalg.norm(v2)
    cos = inner / norms
    rad = np.arccos(np.clip(cos, -1.0, 1.0))
    deg = np.rad2deg(rad)
    return deg

def measure_dihedral_angle(p0,p1,p2,p3):
    """Praxeolitic formula
    1 sqrt, 1 cross product"""

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))

def rmsd(matrix1, matrix2, align=True):
    if align:
        # Flatten matrices and align them based on their centroids
        centroid1 = np.mean(matrix1, axis=0)
        centroid2 = np.mean(matrix2, axis=0)
        matrix1 -= centroid1
        matrix2 -= centroid2

        # Calculate the optimal rotation matrix using Singular Value Decomposition (SVD)
        rotation_matrix, t = kabsch_algorithm(
            matrix2.reshape(3,-1),
            matrix1.reshape(3,-1),
            center=False
        )

        # Apply the rotation to matrix1 and calculate the RMSD
        transformed_matrix1 = (rotation_matrix @ matrix1.T).T
        rmsd_value = np.sqrt(np.mean(np.square(transformed_matrix1 - matrix2)))
    else:
        rmsd_value = np.sqrt(np.mean(np.square(matrix1 - matrix2)))

    return rmsd_value

def rmsd_matrix(matrices):
    num_matrices = len(matrices)
    rmsd_matrix = np.zeros((num_matrices, num_matrices))

    for i in range(num_matrices):
        for j in range(i):  # Only calculate the lower triangular part
            rmsd_matrix[i, j] = rmsd(matrices[i], matrices[j])

    return rmsd_matrix

def rmsd_matrix_compare(matrix_a, matrix_b, align=True):
    rmsd_matrix = np.zeros((len(matrix_a), len(matrix_b)))
    for i in range(len(matrix_a)):
        for j in range(len(matrix_b)):
            rmsd_matrix[i, j] = rmsd(matrix_a[i], matrix_b[j], align)
    return rmsd_matrix

def get_duplicates_rmsd_matrix(matrix, threshold=0.25):
    analysis = np.logical_and(matrix > 0, matrix <= threshold)
    to_delete = np.unique(np.where(analysis)[0])
    return to_delete