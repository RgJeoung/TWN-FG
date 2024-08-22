import os
import sys
import pandas as pd
import numpy as np
import requests
import argparse
from tqdm import tqdm
from bs4 import BeautifulSoup as bs


def DownloadKlifsID(pdb):
    page = requests.get(f'https://klifs.net/api/structures_pdb_list?pdb-codes={pdb}')
    soup = bs(page.text, "html.parser")
    txt = soup.text
    files = txt.split('{')[1:]
    dataframe = pd.DataFrame()
    for file in files:
        lines = file.split(',')
        tmp = {}
        # print(file)
        for line in lines:
            # print(line)
            if line.startswith('"structure_ID"'):
                tmp['ID'] = line[15:]
            elif line.startswith('"kinase"'):
                tmp['kinase'] = line[9:]
            elif line.startswith('"pdb"'):
                tmp['pdb'] = line[6:]
            elif line.startswith('"alt"'):
                tmp['alt'] = line[6:]
            elif line.startswith('"chain"'):
                tmp['chain'] = line[8:]
            elif line.endswith('}') or line.endswith('}]'):
                dataframe = pd.concat([dataframe, pd.DataFrame([tmp])], ignore_index=True)
            else:
                continue
    for idx, ID in enumerate(dataframe['ID']):
        if dataframe['chain'][idx][1:-1] != 'A':
            continue
        else:
            A_chains = dataframe[dataframe['chain'] == 'A']
            if (len(A_chains[A_chains['kinase'] == dataframe['kinase'][idx]]) >= 2) and (
                    dataframe['alt'][idx][1:-1] != 'A'):
                continue
            else:
                return ID


def LoadAbsoluteResidue(pdb):

    # sites = {'AP': [15, 46, 51, 75], 'FP': [10, 51, 72, 81], 'GA': [17, 45, 81], 'SE': [51]} # KinFragLib binding site
    residues = [17, 51]
    ID = DownloadKlifsID(pdb)
    page = requests.get(f'https://klifs.net/details.php?structure_id={ID}')
    soup = bs(page.text, "html.parser")
    lines = str(soup).split('\n')
    real_residues = []
    for line in lines:
        if line.startswith('</tr><tr><td class="residueSearch">'):
            for residue in residues:
                real_res_id = line.find(f'<td class="residueSearch">{residue} ')
                if real_res_id != -1:
                    start = line.find('<span class="xray">', real_res_id + 28) + 19
                    end = line.find('</span>', start)
                    real_residues += [line[start: end]]
        else:
            continue
    return real_residues


def Residue_based_center(inputfile, residue_lst):
    coords = []
    with open(inputfile) as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('ATOM'):
                for residue in residue_lst:
                    if line[22:26].strip() == residue and line[13:15] == 'CA':
                        coord = (float(line[31:38].strip()), float(line[39:46].strip()), float(line[47:54].strip()))
                        coords += [coord]
    center = np.mean(coords, axis=0)
    return tuple(center)


def Klifs_absolute_residue_center(residues, inputfile):
    center = Residue_based_center(inputfile, residues)
    return center


def mkbsfiles(pdb, center, size, inputfile, outputpath):
    site_centers = Klifs_absolute_residue_center(pdb, inputfile)
    # os.mkdir(outputpath + '/' + pdb)

    for site, center in site_centers.items():
        txt_lines = [f'# Boundary file for TWN-Region-Analysis\n# Protein Kinase\n# PDB ID: {pdb}\n\nMETHOD\tpoint\n\nRANGE\t{size}\n\n']
        for i in range(1, 1001):
            txt_lines += [f'COORD\t{str(i)}\t{str(center[0])}\t{str(center[1])}\t{str(center[2])}\n']
        w = open(outputpath + '/' + pdb + '/' + site + '.bd', 'w')
        for txt_line in txt_lines:
            w.write(txt_line)
        w.close()
    return



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Automatic calculation')
    parser.add_argument('-path', '--pdbpath', required=True, help='Set your trajectory directory to make boundary file')
    args = parser.parse_args()

    c_method = input("Select your boundary center method:\n"
                     "[0] Coordinate (x, y, z)\n"
                     "[1] Residue number\n"
                     "[2] Kinase, centre of KLIFS number 17 and 51\n")

    inputpath = args.pdbpath.replace("\\", "/")

    if c_method == '0':  # Give coordinate of center
        center = tuple(float(x) for x in input('Input your coordinate(x, y, z). ex) 1, 2, 3\n').split(', '))
    elif c_method == '1':  # Input residue C alpha chain atom
        residue_number = int(input('Type the residue number to define as a center: ex) 274\n'))
    elif c_method == '2':
        residues = LoadAbsoluteResidue(inputpath.split("/")[-1])
    else:
        sys.exit('Value error: Please write a right c_method.\n')

    size = float(input("Type your boundary size: ex) 20.5\n"))

    txt_lines = [f'# Boundary file for TWN-Region-Analysis\n# PDB ID: {inputpath.split("/")[-1]}\n\nMETHOD  point\n\nRANGE  {size}\n\n']

    for i in tqdm(range(0, len(os.listdir(inputpath + '/a_input')) - 1), desc="Making boundary file..."):
        # define center
        if c_method == '0':  # Give coordinate of center
            pass
        elif c_method == '1':  # Input residue C alpha chain atom
            center = Residue_based_center(inputpath + f'/a_input/{str(i)}.pdb', [residue_number])
        elif c_method == '2':
            center = Klifs_absolute_residue_center(residues, inputpath + f'/a_input/{str(i)}.pdb')
        else:
            sys.exit('Value error: Your given number is not proper or your input file is not appropriate.\n'
                     '             Please check your input file.\n')
        txt_lines += ["COORD  {:6d}  {:8.3f}  {:8.3f}  {:8.3f}\n".format(i, center[0], center[1], center[2])]

    os.mkdir("/".join(inputpath.split("/")[:-2]) + '/boundary')
    os.mkdir("/".join(inputpath.split("/")[:-2]) + '/boundary/' + inputpath.split("/")[-1])

    w = open("/".join(inputpath.split("/")[:-2]) + '/boundary/' + inputpath.split("/")[-1] + '/Center.bd', 'w')
    for txt_line in txt_lines:
        w.write(txt_line)
    w.close()

    print("Boundary file generated.")