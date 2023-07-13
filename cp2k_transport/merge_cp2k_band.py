#!/usr/bin/env python
import sys, linecache
import numpy as np

def merge_kmesh_cp2k(folds):
    n_tot_kpoints = 0
    Kpoints = []
    Eigenvals = []
    
    for fold in folds:
        bs_file = fold + "/cp2k.bs"
        f = open(bs_file,'r')
        firstline = f.readline().strip().split()
        nkpoints, nbands = int(firstline[3]), int(firstline[8])
        n_tot_kpoints += nkpoints
        
        kpoints = np.zeros((nkpoints, 3))
        eigenvals = np.zeros((nbands, nkpoints))
        
        for kpt in range(nkpoints):
            line_num = nkpoints + 2 * (kpt + 1) + nbands * kpt
            kpos = (linecache.getline(bs_file, line_num)).split()[5:8]
            kpoints[kpt][0], kpoints[kpt][1], kpoints[kpt][2] = \
                float(kpos[0]), float(kpos[1]), float(kpos[2])
            Occupation = []
            for band in range(nbands):
                line_num = 2 + nkpoints + 2 * (kpt + 1) + band + nbands * kpt
                ene = float(linecache.getline(bs_file, line_num).split()[1])
                occ = float(linecache.getline(bs_file, line_num).split()[2])
                eigenvals[band][kpt] = ene
                Occupation.append(occ)
    
        Kpoints.append(kpoints)
        Eigenvals.append(eigenvals)
        f.close()
    
    Kpoints = np.vstack(Kpoints)
    Eigenvals = np.transpose(np.hstack(Eigenvals))
    k_weight = 1/n_tot_kpoints
    
    k_point_style = "#  Special point {:<6d}    {:>12.8f}    {:>12.8f}    {:>12.8f}  not specified\n"
    
    band_style = "#  Point {:<6d}    Spin 1:    {:>12.8f}    {:>12.8f}    {:>12.8f}     {:>12.8f}\n#   Band     Energy [eV]      Occupation\n"
    
    eigenval_style = "{:>8d}    {:>12.8f}    {:>12.8f}\n"
    
    
    with open("merged.bs", 'w') as fout:
        fout.write("# Set 1: {} special points, {} k-points, {} bands\n".format(
            n_tot_kpoints, n_tot_kpoints, nbands)
            )
        for k in range(n_tot_kpoints):
            fout.write(
                k_point_style.format(k+1, Kpoints[k][0], Kpoints[k][1], Kpoints[k][2])
            )
    
        for k in range(n_tot_kpoints):
            fout.write(
                band_style.format(k+1, Kpoints[k][0], Kpoints[k][1], Kpoints[k][2], k_weight)
            )
            for band in range(nbands):
                fout.write(
                    eigenval_style.format(band+1, Eigenvals[k][band], Occupation[band])
                    )

if __name__ =="__main__":
    try:
        folds = sys.argv[1:]
        merge_kmesh_cp2k(folds)        
    except ValueError:
        print(f"\nusage: {__file__} foldname*\n")
        sys.exit() 

