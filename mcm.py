from pymaster import NmtWorkspace, NmtField
from os.path import exists, getsize


def mcm_process(mcm_path, mask, b, purify_e, purify_b):
    """ Return MCM workspace for shear - shear
    """
    w = NmtWorkspace()
    if exists(mcm_path) and getsize(mcm_path) > 0:
        w.read_from(mcm_path)
        return w
    f2 = NmtField(mask, [mask, mask], purify_e=purify_e, purify_b=purify_b)
    w.compute_coupling_matrix(f2, f2, b)
    w.write_to(mcm_path)
    return w


def make_pos(mask, ick, b, out_dir):
    """ Return MCM workspace for pos - pos
    """
    path = f"{out_dir}/mcm_gcl_ck{ick}.fits"
    w = NmtWorkspace()
    if exists(path) and getsize(path) > 0:
        w.read_from(path)
        return w
    f1 = NmtField(mask, [mask])
    w.compute_coupling_matrix(f1, f1, b)
    w.write_to(path)
    return w
