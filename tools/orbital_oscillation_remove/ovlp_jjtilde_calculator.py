# import spherical bessel function of the first kind
from scipy.special import spherical_jn as jn
# import simpson integration for rough spherical bessel transform
from scipy.integrate import simps
# import numpy for basic calculation
import numpy as np
# import matplotlib for plotting
import matplotlib.pyplot as plt

#INPUT
bessel_nao_ecut = 100 # in Ry
bessel_nao_rcut = 10 # in a.u.
ecutwfc = 100
#     s  p  d  f
ls = [0, 1, 2, 3] # angular momentum for j_l(qr)
#STRU
# factor from Angstrom to a.u.
lat0 = 1.8897261254578281
cell_dm1 = 10 # dimensionless, real length = cell_dm1 * lat0
cell_dm2 = 10 # dimensionless, real length = cell_dm2 * lat0
cell_dm3 = 10 # dimensionless, real length = cell_dm3 * lat0
latvec = np.array([[cell_dm1, 0, 0], [0, cell_dm2, 0], [0, 0, cell_dm3]])
#ANALYSIS
delta_r_q = 0.01 # to calculate j_l(qr), in a.u.

def generate_truncated_spherical_bessel(ecut, rcut, ls, start_r = 0, delta_r = 0.0001, delta_r_q = 0.01, plot_jntilde = False):
    """
    Generate truncated spherical Bessel functions with given cutoff radius, angular momentums
    """
    nq = int(np.sqrt(ecut) * rcut / np.pi)
    qs = [] # qs for every angluar momentum
    for l in ls:
        roots = []
        r = start_r
        while len(roots) < nq:
            if jn(l, r)*jn(l, r+delta_r) < 0:
                roots.append(r)
            r += delta_r
        #print("l = ", l)
        #print("roots = ", roots)

        # then we suppose for every q in qs, j_l(qrcut) = 0, it means q[i]*rcut = roots[i], q[i] = roots[i]/rcut
        r = np.arange(0, rcut, delta_r_q)
        qs_l = np.array(roots) / rcut

        if plot_jntilde:
            for q in qs_l:
                jlq = jn(l, q*r)
                plt.plot(r, jlq, label="l = "+str(l)+", q = "+str(round(q, 4)))
            plt.xlabel("r (a.u.)", fontsize=14)
            plt.ylabel("j_l(q*r)", fontsize=14)
            plt.title("generated truncated spherical bessel function", fontsize=14)
            plt.legend()
            plt.show()

        qs.append(qs_l)
    
    return qs

def generate_g_vector(lat0, latvec, ecutwfc):
    GT = np.linalg.inv(latvec) # it is not a good idea to use all capital letters for variable name
    G = np.transpose(GT)
    GGT = np.dot(G, GT)

    tpiba = 2 * np.pi / lat0 # in a.u.-1

    radius_kspace = np.sqrt(ecutwfc) # ekin = 0.5 * k**2 = 0.5 * ecutwfc
                                                # therefore k = sqrt(ecutwfc), k in a.u.-1, ecutwfc in Ry
    arg = [
        int(radius_kspace/tpiba*np.linalg.norm(latvec[0])),
        int(radius_kspace/tpiba*np.linalg.norm(latvec[1])),
        int(radius_kspace/tpiba*np.linalg.norm(latvec[2]))
    ]
    g_direct_vectors = []
    g_cartesian_vectors = []
    #kpt_max = np.array([ibox[0]+1, ibox[1]+1, ibox[2]+1])
    #print("kpt_max = ", kpt_max)
    for i in range(-arg[0], arg[0]+1):
        for j in range(-arg[1], arg[1]+1):
            for k in range(-arg[2], arg[2]+1):
                kpt = np.array([i, j, k])
                if np.dot(kpt, np.dot(GGT, kpt))*(tpiba**2) <= ecutwfc:
                    g_direct_vectors.append(kpt)
                    g_cartesian_vectors.append(np.dot(G, kpt) * tpiba)
                else:
                    #print("kpt = ", kpt, " is not included in the k-space")
                    pass
    g_direct_vectors = np.array(g_direct_vectors)
    g_cartesian_vectors = np.array(g_cartesian_vectors)

    return g_direct_vectors, g_cartesian_vectors

if __name__ == "__main__":
    
    """ generate g-vectors and truncated spherical bessel functions """

    g_direct_vectors, g_cartesian_vectors = generate_g_vector(lat0, latvec, ecutwfc)
    print("number of g-vectors: ", len(g_direct_vectors))
    qs = generate_truncated_spherical_bessel(bessel_nao_ecut, bessel_nao_rcut, ls, plot_jntilde=True)
    #g_norm = np.linalg.norm(g_cartesian_vectors[-1])
    g_norms = []
    g_norms_interval = 100
    for ignorm in range(len(g_cartesian_vectors)):
        if ignorm%g_norms_interval == 0:
            g_norms.append(np.linalg.norm(g_cartesian_vectors[ignorm]))

    #print("at |g| = ", g_norm, " a.u.-1")
    # generate grid points
    r = np.arange(0, bessel_nao_rcut, delta_r_q)

    # enable sorting of g_norms
    g_norms = np.sort(g_norms)
    jlgs = []
    # generate j_l(g*r) in spherical expansion of plane wave of g
    for g_norm in g_norms:
        jlgs_g = []
        for l in ls:
            jlg = jn(l, g_norm*r)
            jlgs_g.append(jlg)
        jlgs.append(jlgs_g)
    print("will analyze g-vectors at gamma point, g-vectors' norm are: (in a.u.-1)")
    print(g_norms)
    print("number of g-vectors to analyze: ", len(g_norms))
    # calculate the overlap between j_l(g*r) and j_l(q*r) in range 0 to rcut

    """ calculate overlap between j_l(g*r) and j_l(q*r) for every q, l, g """

    norm_ovlp = []
    ovlp = []
    for il, l in enumerate(ls):
        norm_ovlp_l = []
        ovlp_l = []
        print("for angular momentum l = ", l, ", will analyze q-vectors: (in a.u.-1)")
        for q in qs[il]:
            ovlps_lq = []
            #sum_ovlp = []
            for ig in range(len(g_norms)):
                g_norm = g_norms[ig]
                jqr = jn(l, q*r)
                ovlp_lqg = simps(jqr*jlgs[ig][il]*r**2, r)
                ovlps_lq.append(ovlp_lqg)
            # draw overlap between j_l(g*r) and j_l(q*r) for current q, l, variable is |G|
            plt.plot(g_norms, np.log10(np.abs(ovlps_lq)), '-', label="l = "+str(l)+", |q| = "+str(round(q, 4))+" a.u.-1")
            norm_ovlp_l.append(np.linalg.norm(ovlps_lq))
            ovlp_l.append(ovlps_lq)

            print("sum of overlaps of present q over all g-vectors: q = ", round(q, 4), ", sum = ", round(norm_ovlp_l[-1], 4))
        norm_ovlp.append(norm_ovlp_l)
        ovlp.append(ovlp_l)
    
    # add y = 0 line
    plt.plot([0, g_norms[-1]], [0, 0], 'k--')
    plt.xlabel("G (a.u.-1)", fontsize=20)
    plt.ylabel(r'$log_\{10\}[\langle j_l(g)|\tilde{j}_l(q)\rangle(l, q; G)]$', fontsize=20)
    # add legend in two columns
    plt.legend(ncol = 3, fontsize = 14)
    plt.show()

    """ plot sum of overlaps over all g-vectors """

    for il, l in enumerate(ls):
        plt.plot(qs[il], norm_ovlp[il], '-', label="l = "+str(l), linewidth = 3)
    # again, add y = 0 line
    plt.plot([0, qs[0][-1]], [0, 0], 'k--')
    plt.xlabel("q (a.u.-1)", fontsize=20)
    plt.ylabel("norm of overlaps over all g-vectors", fontsize=20)
    plt.legend(fontsize = 14)
    plt.show()

    """ plot hotmap of overlaps per l, overlap dimension (nq, nG)"""
    """
    # plot hotmap of logrithm of abs of overlaps per l, overlap dimension (nq, nG)
    # draw each l's as one subplot
    _log10_abs_ovlp = []
    for il in range(len(ls)):
        for iq in range(len(qs[il])):
            for ig in range(len(g_norms)):
                _log10_abs_ovlp.append(np.log10(np.abs(ovlp[il][iq][ig])))
    _log10_abs_ovlp = np.array(_log10_abs_ovlp)
    ovlp = np.reshape(_log10_abs_ovlp, (len(ls), len(qs[0]), len(g_norms)))
    # plot
    ncol = 2
    nrow = int(len(ls)/ncol)
    fig, axs = plt.subplots(nrow, ncol, sharex=True, sharey=True, squeeze=False)
    # set size
    fig.set_size_inches(10, 10)
    fig.suptitle('log10 of abs of overlaps per l', fontsize=20)
    irow = 0
    icol = 0
    for il, l in enumerate(ls):
        print("dimension of ovlp[", il, "] = ", ovlp[il].shape)
        # set the same color scale for all subplots
        axs[irow][icol].imshow(ovlp[il], 
                               cmap='plasma', 
                               interpolation='nearest', 
                               origin='lower', 
                               extent=[g_norms[0], g_norms[-1], qs[il][0], qs[il][-1]],
                               aspect='auto')
        
        axs[irow][icol].set_xlabel("G (a.u.-1)", fontsize=20)
        axs[irow][icol].set_ylabel("q (a.u.-1)", fontsize=20)
        axs[irow][icol].set_title('l = '+str(l))
        icol += 1
        if icol == ncol:
            icol = 0
            irow += 1
    plt.show()
    """

    """plot hotmap of overlaps per l, overlap dimension (nq, nG), use contourf"""

    # plot hotmap of logrithm of abs of overlaps per l, overlap dimension (nq, nG)
    # draw each l's as one subplot
    _log10_abs_ovlp = []
    for il in range(len(ls)):
        for iq in range(len(qs[il])):
            for ig in range(len(g_norms)):
                _log10_abs_ovlp.append(np.log10(np.abs(ovlp[il][iq][ig])))
    _log10_abs_ovlp = np.array(_log10_abs_ovlp)
    ovlp = np.reshape(_log10_abs_ovlp, (len(ls), len(qs[0]), len(g_norms)))
    # plot
    ncol = 2
    nrow = int(len(ls)/ncol)
    fig, axs = plt.subplots(nrow, ncol, sharex=True, sharey=True, squeeze=False)
    # set size
    fig.set_size_inches(10, 10)
    fig.suptitle('log10 of abs of overlaps per l', fontsize=20)
    irow = 0
    icol = 0
    # set x axis as g_norms, y axis as qs
    x = g_norms
    y = []
    for il in range(len(ls)):
        y.append(qs[il])
    for il, l in enumerate(ls):
        print("dimension of ovlp[", il, "] = ", ovlp[il].shape)
        # set the same color scale for all subplots
        vmax = np.max(ovlp[il])
        vmin = -10
        print("contourf: l = ", l, "vmax = ", vmax, ", vmin = ", vmin)
        axs[irow][icol].contourf(x, y[il], ovlp[il], 
                               cmap='plasma', 
                               interpolation='nearest', 
                               origin='lower', 
                               extent=[x[0], x[-1], y[il][0], y[il][-1]],
                               aspect='auto',
                               vmin=vmin,
                               vmax=vmax)
        
        axs[irow][icol].set_xlabel("G (a.u.-1)", fontsize=20)
        axs[irow][icol].set_ylabel("q (a.u.-1)", fontsize=20)
        axs[irow][icol].set_title('l = '+str(l))
        icol += 1
        if icol == ncol:
            icol = 0
            irow += 1
    plt.show()