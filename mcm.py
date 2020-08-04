from pymaster import NmtWorkspace, NmtField
from os.path import exists


def mcm_process(mcm_path, mask, b, purify_e, purify_b):
    w = NmtWorkspace()
    if exists(mcm_path):
        w.read_from(mcm_path)
        return w
    f2 = NmtField(mask, [mask, mask], purify_e=purify_e, purify_b=purify_b)
    w.compute_coupling_matrix(f2, f2, b)
    w.write_to(mcm_path)
    return w
