import numpy as np


def lazy_property(fn):
    """
    Implementation borrowed from https://towardsdatascience.com/what-is-lazy-evaluation-in-python-9efb1d3bfed0
    """
    attr_name = "_lazy_" + fn.__name__

    @property
    def _lazy_property(self):
        if not hasattr(self, attr_name):
            setattr(self, attr_name, fn(self))
        return getattr(self, attr_name)

    return _lazy_property


def n_to_subtract(atoms_removed_inds, atom_ind):
    return int(np.argwhere(np.sort(atoms_removed_inds + [atom_ind]) == atom_ind)[0][0])
