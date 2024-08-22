import os
import logging
import argparse
import numpy as np
from math import dist
from pathlib import Path


# module info option
logger = logging.getLogger(__name__)


# set logger
def set_log(path_output, log_message):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(path_output / log_message),
            logging.StreamHandler()
        ]
    )


def sdf_reader(sdf_file):
    point_dicts = {}
    line_dicts = {}
    with open(sdf_file, 'r') as f:
        lines = f.readlines()
        naming = True
        points = []
        line_for_1pattern = []
        for line in lines:
            line_for_1pattern += [line]
            if naming:
                name = line.strip()
                naming = False
            if len(line) == 40 and len(line.split()) == 6 and line.split()[3] != 'H' and line.split()[3] != '?':
                points += [(float(line.split()[0]), float(line.split()[1]), float(line.split()[2]))]
            if line.startswith('$$$$'):
                point_dicts[name] = points
                points = []
                line_dicts[name] = line_for_1pattern
                line_for_1pattern = []
                naming = True
    return point_dicts, line_dicts


# trans format
def trans_format(form, *vals):
    # sdf format
    sdf_header = "{:s}\n  TWNREGION          3D                             0\n"
    sdf_count = "{:3d}{:3d}  0  0  0  0  0  0  0  0999 V2000"
    sdf_atom = "{:10.4f}{:10.4f}{:10.4f} {:1s} {:3d}{:3d}"
    sdf_bond = "{:3d}{:3d}{:3d}  0  0  0"
    sdf_prop_s = "> <{:s}>\n{:s}\n"
    sdf_prop_d = "> <{:s}>\n{:d}\n"
    #pdb_format
    pdb_line = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}"
    pdb_ter = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}"
    pdb_model = '{:6s}     {:>3d}'
    trans = {'sdf_header': sdf_header, 'sdf_count': sdf_count, 'sdf_atom': sdf_atom,
             'sdf_bond': sdf_bond, 'sdf_prop_s': sdf_prop_s, 'sdf_prop_d': sdf_prop_d,
             'pdb_line': pdb_line, 'pdb_ter': pdb_ter, 'pdb_model': pdb_model}
    return trans[form].format(*vals)


def region_extractor(inputpath, protein):
    os.mkdir(inputpath + f'/TWN-Region/{protein}')
    file = inputpath + f'/TWN-Pattern/{protein}/TWN.sdf'
    point_dicts, line_dicts = sdf_reader(file)
    ring_type = int(list(line_dicts.values())[0][3].split()[0])
    pattern_counting = {}
    for name, points in point_dicts.items():
        mean_points = np.mean(points, axis=0)
        pattern_counting[name] = []
        for n, p in point_dicts.items():
            if dist(mean_points, np.mean(p, axis=0)) <= 1:
                pattern_counting[name] += [n]

    sorted_pattern_counting = dict(sorted(pattern_counting.items(), key=lambda x: len(set(sum([line_dicts[x][ring_type + 9].strip().split() for x in x[1]], []))), reverse=True))
    region_number = 1
    registered_patterns = []
    for g_name, g_patterns in sorted_pattern_counting.items():
        unique = True
        for g_pattern in g_patterns:
            if g_pattern in registered_patterns:
                unique = False
        if unique:
            region_name = f'TWN_Region_{region_number}'
            frequency = len(set(sum([line_dicts[x][ring_type + 9].strip().split() for x in g_patterns], [])))
            w= open(inputpath + f'/TWN-Region/{protein}/{region_name}.sdf', 'w')
            header_row = trans_format('sdf_header', 'TWN_Region_' + str(region_number))

            count_row = trans_format('sdf_count', ring_type * len(g_patterns), 0)
            w.write(header_row + '\n' + count_row + '\n')
            for g_pattern in g_patterns:
                for g_line in line_dicts[g_pattern][4:4+ring_type]:
                    w.write(g_line)
            w.write(f'M  END\n\n> <region.frequency>\n{frequency}\n\n> <twn.pattern.names>\n{"  ".join(g_patterns)}\n')
            w.close()
            region_number += 1
            registered_patterns += g_patterns
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Identify TWN Regions by its union frequency.')
    parser.add_argument('-d', '--directory', required=True, help='Set your trajectory directory')
    args = parser.parse_args()

    inputpath = args.directory.replace("\\", "/")
    proteins = os.listdir(inputpath + "/TWN-Pattern")
    os.mkdir(inputpath + "/TWN-Region")
    for protein in proteins:
        log_path = Path(inputpath + '/logs/' + protein)
        set_log(log_path, "TWN-Region-Analysis.log")
        logger.info(f"Start Region identification...")
        region_extractor(inputpath, protein)
        logger.info(f"Region identification complete.")
        logger.info(f"{len(os.listdir(inputpath + '/TWN-Region/' + protein))} Regions are extracted from {protein}.")
    logger.info(f"Process finished.")