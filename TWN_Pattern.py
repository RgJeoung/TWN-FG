import logging
import argparse
import numpy as np
from tqdm import tqdm
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


# trans format
def trans_format(form, *vals):
    # sdf format
    sdf_header = "{:s}\n  TWNPATTERN          3D                             0\n"
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


def boundary_reader(boundary_file):
    boundary_inform = []
    with open(boundary_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('METHOD'):
                boundary_inform += ["_".join(line.split()[1:]).rstrip()]
            if line.startswith('RANGE'):
                boundary_inform += [float(line.split()[1].strip())]
            if line.startswith('COORD'):
                boundary_inform += [[int(line.split()[1]), tuple((float(line.split()[2].strip()), float(line.split()[3].strip()), float(line.split()[4].strip())))]]
            if line.startswith('RESIDUE'):
                boundary_inform += ["-".join(line.split()[1:]).rstrip()]
    return boundary_inform


def pdb_spliter(pdb_line):
    type = pdb_line[0:6].strip()
    atm_index = pdb_line[6:11].strip()
    atm_type = pdb_line[11:17].strip()
    residue_name = pdb_line[17:20].strip()
    chain = pdb_line[21:22].strip()
    residue_index = pdb_line[22:26].strip()
    x_coord = pdb_line[30:38].strip()
    y_coord = pdb_line[38:46].strip()
    z_coord = pdb_line[46:54].strip()
    frequency = pdb_line[54:60].strip()
    b_factor = pdb_line[60:66].strip()
    element = pdb_line[66:].strip()
    return [type, atm_index, atm_type, residue_name, chain, residue_index,
            x_coord, y_coord, z_coord, frequency, b_factor, element]


# read protein
def trajectory_reader(trajectory_file, boundary_file):
    # boundary setting
    boundary_inform = boundary_reader(boundary_file)
    centers = []
    residue_centers = []
    if boundary_inform[0] != 'residue_center_extraction':
        for trj, coord in boundary_inform[2:]:
            if trajectory_file.stem == str(trj):
                centers += [coord]
    # get protein inform
    with open(trajectory_file, 'r') as f:
        lines = f.readlines()
        protein_inform = []
        for line in lines:
            line = line.rstrip()
            atoms = pdb_spliter(line)

            if atoms[0] != 'ATOM':
                protein_inform += [line]
                continue

            elif atoms[3] == 'SOL' or atoms[3] == 'CL' or atoms[3] == 'NA':
                continue

            else:
                if boundary_inform[0] == 'residue_center_extraction' and (atoms[3] + atoms[5]) in boundary_inform[2].split("-") and atoms[2] == 'CA':
                    residue_centers += [tuple((float(atoms[6]), float(atoms[7]), float(atoms[8])))]
                pdb_row = trans_format('pdb_line', atoms[0], int(atoms[1]), atoms[2], '', atoms[3], atoms[4], int(atoms[5]), '',
                                       float(atoms[6]), float(atoms[7]), float(atoms[8]),
                                       float(atoms[9]), float(atoms[10]), atoms[11], '')
                protein_inform += [pdb_row]
        if len(residue_centers) != 0:
            centers = [np.mean(residue_centers, axis=0)]
        water_inform = []
        first_atm = True
        ht = 0
        for line in lines:
            line = line.rstrip()
            atoms = pdb_spliter(line)

            if atoms[3] == 'SOL':
                # save atom number
                if int(atoms[1]) == 0:
                    ht += 100000
                atm_num = int(atoms[1]) + ht
                # define first atom
                if first_atm:
                    first_num = atm_num
                    first_atm = False
                if atoms[2] == 'OW':
                    w_num = int((atm_num-first_num)/3 + 1)
                    if any([dist((float(atoms[6]), float(atoms[7]), float(atoms[8])), center) <= boundary_inform[1] for center in centers]):
                        water_row = trans_format('pdb_line', atoms[0], int(atoms[1]), atoms[2], '', 'W' + "{:0>6d}".format(w_num)[:2], '', int("{:0>6d}".format(w_num)[2:]), '',
                                                 float(atoms[6]), float(atoms[7]), float(atoms[8]),
                                                 float(atoms[9]), float(atoms[10]), atoms[11], '')
                        water_inform += [water_row]
    return protein_inform, water_inform


# read single water
def single_reader(trj, boundary_file):
    single_box = {}
    for trajectory_file in tqdm(trj):
        trajectory = trajectory_file.stem
        file_box = []
        lines = trajectory_reader(trajectory_file, boundary_file)[1]
        for line in lines:
            atoms = pdb_spliter(line)
            if atoms[2] != 'OW':
                continue
            if len(atoms) == 11:
                atoms += ['O']
            if len(atoms) != 12:
                continue
            single_data = [atoms[1], atoms[2], trajectory, atoms[3] + atoms[5],
                        float(atoms[6]), float(atoms[7]), float(atoms[8]),
                        float(atoms[9]), float(atoms[10]), atoms[11]]
            file_box += [single_data]
        single_box[trajectory] = file_box
    return single_box


# read twn water
def TWN_reader(twn):
    twn_box = {}
    twn_atom_num = 1
    for twn_file in tqdm(twn):
        twn_name = twn_file.stem
        file_box = []
        with open(twn_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.rstrip()
                atoms = pdb_spliter(line)
                if atoms[2] != 'OW':
                    continue
                if len(atoms) == 11:
                    atoms += ['O']
                if len(atoms) != 12:
                    continue
                twn_data = [twn_atom_num, atoms[2], twn_name, atoms[3] + atoms[5],
                            float(atoms[6]), float(atoms[7]), float(atoms[8]),
                            float(atoms[9]), float(atoms[10]), atoms[11]]
                file_box += [twn_data]
            twn_atom_num += 1
            twn_box[twn_name] = file_box
    return twn_box


# write twn water
def TWN_writer(twn_box, single_box, twn_out):
    logger.info(f'Start identifying TWN-Patterns...')
    TWN_patterns = {}
    for twn_name, twn_data in twn_box.items():
        tmp_TWN_pattern = {}
        TWN_patterns[twn_name] = {}
        for cluster_center in twn_data:
            tmp_TWN_pattern[cluster_center[3]] = {}
            for trj, single_data in single_box.items():
                if trj != twn_name.split("_")[1]:
                    for single_water in single_data:
                        if dist(tuple(cluster_center[4:7]), tuple(single_water[4:7])) <= 1.0:
                            tmp_TWN_pattern[cluster_center[3]][single_water[2]] = single_water[3]
        accepted_trj = {key: True for key in list(set([x for x in single_box.keys()]))}
        accepted_trj_copy = list(accepted_trj.keys())
        for tmp_trj in tmp_TWN_pattern.values():
            for trajec in accepted_trj_copy:
                if trajec not in list(tmp_trj.keys()):
                    accepted_trj[trajec] = False
        TWN_patterns[twn_name][twn_name.split("_")[1]] = "-".join(list(tmp_TWN_pattern.keys()))
        for atk, atv in accepted_trj.items():
            if atv:
                TWN_patterns[twn_name][atk] = "-".join([v[atk] for k, v in tmp_TWN_pattern.items()])

    sorted_TWN_patterns = dict(sorted(TWN_patterns.items(), key=lambda x: len(list(x[1].keys())), reverse=True))
    logger.info(f'Total TWN is {len(sorted_TWN_patterns.keys())}.')
    logger.info(f'Start extracting unique region frequent TWN Patterns...')
    sd_idx = 0
    first_twn = True
    pattern_box = {}
    twn_each = []
    for t_name, t_trjs in tqdm(sorted_TWN_patterns.items()):
        unique = False
        if first_twn:
            unique = True
        else:
            tmp_box = pattern_box
            unique_check = True
            for twn_n, twn_c in tmp_box.items():
                unique_pass = True
                for twn_cn in twn_c:
                    unique_water = True
                    twn_v = [x for x in twn_box[twn_n] if x[3] == twn_cn]
                    twn_coord = twn_v[0][4:7]
                    for twn_d in twn_box[t_name]:
                        if dist(twn_coord, tuple(twn_d[4:7])) <= 1.0:
                            unique_water = False
                            break
                    if unique_water:
                        unique_pass = False
                if unique_pass:
                    unique_check = False
            if unique_check:
                unique = True
        if unique:
            header_row = trans_format('sdf_header', 'TWN_Pattern_' + str(sd_idx + 1))
            sd_idx += 1
            twn_props = {'twn.center.name': t_name, 'twn.occupation.trjs': "  ".join(list(t_trjs.keys())),
                         'twn.frequency': len(list(t_trjs.keys())),
                         'twn.w.names': "  ".join(list(t_trjs.values()))}
            count_row = trans_format('sdf_count', len(twn_box[t_name]), 0)
            if twn_props['twn.frequency'] >= 2:
                twn_each += [[header_row] + [count_row] +
                             [trans_format('sdf_atom', x[4], x[5], x[6], "O", 0, 0) for x in list(twn_box[t_name])] +
                             ['M  END'] + [trans_format('sdf_prop_s', 'twn.center.name', twn_props['twn.center.name'])] +
                             [trans_format('sdf_prop_s', 'twn.occupation.trjs', twn_props['twn.occupation.trjs'])] +
                             [trans_format('sdf_prop_d', 'twn.frequency', twn_props['twn.frequency'])] +
                             [trans_format('sdf_prop_s', 'twn.w.names', twn_props['twn.w.names'])] + ['$$$$']]
                pattern_box[t_name] = twn_props["twn.w.names"].split("  ")[0].split("-")
            else:
                break
            first_twn = False
    logger.info(f'Unique TWN Patterns = {len(sorted_TWN_patterns.keys())}(Total TWN Patterns) - {len(sorted_TWN_patterns.keys()) - len(pattern_box.keys())}(Duplicated TWN Patterns) = {len(pattern_box.keys())}')
    twn_sdf = sum(twn_each, [])
    # write twn
    with open(twn_out, 'w') as f:
        for line in twn_sdf:
            f.write(f"{line}\n")

    return twn_sdf


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Identify TWN Patterns by its frequency.')
    parser.add_argument('-trj', '--trajectory', required=True, help='Set your trajectory directory')
    parser.add_argument('-bd', '--boundary', required=True, help='Set your boundary file')
    parser.add_argument('-twn', '--twn_water', required=True, help='Set your twn directory')
    parser.add_argument('-o', '--output', required=True, help='Set your output directory')
    parser.add_argument('-l', '--log', required=True, help='Set your log directory')
    args = parser.parse_args()

    # set output
    output_path = Path(rf"{args.output}")
    log_path = Path(rf"{args.log}")
    set_log(log_path, "TWN-Region-Analysis.log")

    logger.info(f"Analysis for {log_path.stem}")
    logger.info(f'Start Pattern identification...')
    # set twn
    logger.info(f'Loading TWN data...')
    TWN_path = Path(rf"{args.twn_water}")
    TWN = [twn for twn in TWN_path.glob('./*.pdb')]
    TWN_out = output_path / f"TWN.sdf"
    TWN_box = TWN_reader(TWN)
    logger.info(f'Set TWN water : {TWN_path} | {len(TWN)} files in folder')

    # set single water
    logger.info(f'Loading single water data...')
    single_path = Path(rf"{args.trajectory}")
    single_water = [single for single in single_path.glob('./*.pdb')]
    single_box = single_reader(single_water, args.boundary)
    logger.info(f'Set single water : {single_path} | {len(single_water)} files in folder')
    logger.info(rf'Limited boundary: Use {args.boundary}')

    # TWN pattern identification
    TWN_writer(TWN_box, single_box, TWN_out)
    logger.info(f'Saved pdb : {TWN_out}')

    logger.info('Pattern identification complete.')