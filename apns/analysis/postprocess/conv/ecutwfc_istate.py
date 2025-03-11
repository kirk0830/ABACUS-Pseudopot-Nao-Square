import apns.analysis.postprocess.read_abacus_out as amarao
import apns.analysis.postprocess.conv.kernel as amack
import numpy as np
import os

def search_istate(search_domain: str):
    # istate is indiced as [ispin][ik][iband][icol],
    # in which the first col is index of band starting from 1
    # the second col is eigenvalue in Ry
    # the third col is occupation number
    path, istates = [], []
    for root, dirs, files in os.walk(search_domain):
        if "istate.info" in files:
            # get the istate
            istate, _ = amarao.read_istate(os.path.join(root, "istate.info"))
            if istate is not None:
                path.append(root)
                istates.append(istate)
    return path, istates

def cal_band_diff(istate, istate_ref):
    """istate specific calculator, to pass as argument to amack.calculate
    return one scalar value.
    
    istate: band structure at one ecutwfc
    istate_ref: band structure at the reference ecutwfc
    
    Formula:
    $\eta_v(A, B) = min_{\omega}*\sqrt{
        \frac{\sum_{nk}{
            \tilde{f}_{nk}(e^A_{nk}-e^B_{nk}+\omega)^2}
            }{\sum_{nk}{\tilde{f}_{nk}}}
    }$

    Here for one pseudopotential, it is simplified that omega is not
    needed to be optimized, so the formula is simplified as:
    $\eta_v(A, B) = \sqrt{
        \frac{\sum_{nk}{
            \tilde{f}_{nk}(e^A_{nk}-e^B_{nk})^2}
            }{\sum_{nk}{\tilde{f}_{nk}}}
    }$

    istate and istate_ref are in the form of [ispin][ik][iband][icol],
    f_{nk} therefore should sum across spin, k, and band.
    """
    val = 0.0 # the value to be returned

    # not sure if istate and istate_ref is already np.array. 
    istate = np.array(istate)
    istate_ref = np.array(istate_ref)

    mode = "nonsence"
    if len(istate.shape) - len(istate_ref.shape) == 1:
        # this means istate is a list in which each element should have
        # the same shape as istate_ref, therefore assert
        assert istate.shape[1:] == istate_ref.shape
        mode = "vector"
    elif len(istate.shape) - len(istate_ref.shape) == 0:
        assert istate.shape == istate_ref.shape
        mode = "scalar"
    else:
        mode = "nonsence"
        raise ValueError("Unexpected value to evaluate istate convergence.")

    assert mode != "nonsence"
    istate = np.array([istate]) if mode == "scalar" else istate
    vals = []
    for i in range(istate.shape[0]):
        f_nk = istate[i, :, :, :, 2].flatten() # [ispin*ik*iband]
        f_nk_ref = istate_ref[:, :, :, 2].flatten() # [ispin*ik*iband]
        # element-wise geometrical average on f_nk and f_nk_ref
        f_nk = np.sqrt(f_nk * f_nk_ref)

        e_nk = istate[i, :, :, :, 1].flatten() # [ispin*ik*iband]
        e_nk_ref = istate_ref[:, :, :, 1].flatten() # [ispin*ik*iband]

        vals.append(np.sqrt(np.sum(f_nk * (e_nk - e_nk_ref)**2) / np.sum(f_nk)))

    return vals if mode == "vector" else vals[0]

def run(search_domain: str, thr: float = 1e-2):
    first, second = amack.search(search_domain=search_domain,
                                 searcher=search_istate,
                                 scalarizer=cal_band_diff)
    conv = amack.calculate(first, second, thr=thr)
    return conv, second

import unittest
class TestEcutwfcIstate(unittest.TestCase):

    def test_calculator(self):
        # nspin = 1
        istate_1 = [[
           # ib, e,   f
            [[1, 0.0, 1.0], 
             [2, 1.0, 1.0]], # ik = 0
            [[1, 0.0, 1.0], 
             [2, 1.0, 1.0]]  # ik = 1
        ]]
        istate_2 = [[
           # ib, e,   f
            [[1, 0.0, 1.0], 
             [2, 1.0, 1.0]], # ik = 0
            [[1, 0.0, 1.0], 
             [2, 1.0, 1.0]]  # ik = 1
        ]]
        istate_1 = np.array(istate_1)
        istate_2 = np.array(istate_2)
        val = cal_band_diff(istate_1, istate_2)
        self.assertAlmostEqual(val, 0.0)
        # not identical istate
        istate_2 = [[
              # ib, e,   f
                [[1, 0.0, 1.0], 
                 [2, 0.5, 1.0]], # ik = 0
                [[1, 0.0, 1.0], 
                 [2, 0.3, 1.0]]  # ik = 1
          ]]
        istate_2 = np.array(istate_2)
        val = cal_band_diff(istate_1, istate_2)
        self.assertGreater(val, 0.0)

        # nspin = 2
        istate_1 = [[
           # ib, e,   f
            [[1, 0.0, 1.0], 
             [2, 1.0, 1.0]], # ik = 0
            [[1, 0.0, 1.0], 
             [2, 1.0, 1.0]]  # ik = 1
        ],
        [
           # ib, e,   f
            [[1, 0.0, 1.0], 
             [2, 1.0, 1.0]], # ik = 0
            [[1, 0.0, 1.0], 
             [2, 1.0, 1.0]]  # ik = 1
        ]]
        istate_2 = [[
           # ib, e,   f
            [[1, 0.0, 1.0], 
             [2, 1.0, 1.0]], # ik = 0
            [[1, 0.0, 1.0], 
             [2, 1.0, 1.0]]  # ik = 1
        ],
        [
           # ib, e,   f
            [[1, 0.0, 1.0], 
             [2, 1.0, 1.0]], # ik = 0
            [[1, 0.0, 1.0], 
             [2, 1.0, 1.0]]  # ik = 1
        ]]
        istate_1 = np.array(istate_1)
        istate_2 = np.array(istate_2)
        val = cal_band_diff(istate_1, istate_2)
        self.assertAlmostEqual(val, 0.0)
        # not identical istate
        istate_2 = [[
              # ib, e,   f
                [[1, 0.0, 1.0], 
                 [2, 0.4, 1.0]], # ik = 0
                [[1, 0.1, 1.0], 
                 [2, 0.5, 1.0]]  # ik = 1
          ],
          [
              # ib, e,   f
                [[1, 0.0, 1.0], 
                 [2, 0.5, 1.0]], # ik = 0
                [[1, -0.2, 1.0], 
                 [2, 0.7, 1.0]]  # ik = 1
          ]]
        istate_2 = np.array(istate_2)
        val = cal_band_diff(istate_1, istate_2)
        self.assertGreater(val, 0.0)

if __name__ == "__main__":

    path = "../11668634"
    result = run(path)
    print(result)
    exit()
    unittest.main()