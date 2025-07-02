#!/usr/bin/env python

import os,sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import lines

def read_elevels(filename):
    pwd = os.getcwd()
    filename = os.path.join(pwd,filename)
    if os.path.isfile(filename):
        os.system(f'grep "MO|" {filename}> mo_eigenvals.txt')
    else:
        print(f'{filename} not find!')
        os.sys.exit(0)
    f = open('mo_eigenvals.txt', 'r')

    lines = f.readlines()
    targetline = 0
    for linenum, line in enumerate(lines):
        line_str = lines[linenum].strip().split()
        if 'EIGENVALUES' in line_str: 
            targetline = linenum
            
    eigenval = []
    occ = []
    for linenum in range(targetline+3, len(lines)):
       line_str = lines[linenum].strip().split()
       if 'Sum:' in line_str:
           efermi = lines[linenum+1].strip().split()[4]
           return eigenval, float(efermi), occ
       eigenval.append(float(line_str[3]))
       occ.append(float(line_str[4]))

def plot_elevels(e_levels,efermi, erange, outfile='elevels.png'):
    e_levels = np.array(e_levels)
    
    padding = (erange[1] - erange[0]) / 50 
    figure = plt.figure(figsize=(2,3))
    ax = plt.gca()
    ax.get_xaxis().set_visible(False) 
    ax.get_yaxis().set_label_text('Energy')
    ax.set_ylim(erange[0] - padding, erange[1] + padding) 
    for eval_ in e_levels:
        if erange[0] < eval_ < efermi:
            line = lines.Line2D((0.3, 0.7), (eval_, eval_), c='darkblue') 
            ax.add_line(line)
        if efermi < eval_ < erange[1]:
            line = lines.Line2D((0.3, 0.7), (eval_, eval_), c='orange') 
            ax.add_line(line)
    plt.tight_layout()
    plt.savefig(outfile,dpi=300)

def _main():
    if len(sys.argv) < 2:
        print("    usage: plot_cp2k_elevels.py outfilename emin emax [eref]\n")
        os.sys.exit()
    eigenval,efermi,occ = read_elevels(sys.argv[1])        
    erange = (float(sys.argv[2]), float(sys.argv[3]))
    try: 
        eref = float(sys.argv[4])
    except IndexError:
        eref = 0
        pass
    efermi -= eref
    eigenval = [ e - eref  for e in eigenval ]
    plot_elevels(eigenval,efermi,erange)

######################
_main()
