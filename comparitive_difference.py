class comp_diff(object):
    """The class comp_diff is the normalized comparison difference，
    analytical_phi is the analytical potential of an object,
    multipole_phi is the multipole expansion approximation potential."""

    def __init__(self, analytical_phi, multipole_phi):
        self.difference = (multipole_phi-analytical_phi)/analytical_phi
