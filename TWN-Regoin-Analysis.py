import os
import sys
import argparse
import subprocess
from pathlib import Path


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Automatic calculation')
    parser.add_argument('-d', '--directory', required=True, help='Set your directory to analyze')
    args = parser.parse_args()

    # Pattern identification
    path = f"{args.directory}".replace("\\", "/")
    proteins = os.listdir(path + "/TWN")
    trajectory_path = path + "/trajectory"
    boundary_path = path + "/boundary"
    output_path = path + "/TWN-Pattern"
    log_path = path + "/logs"
    os.mkdir(output_path)
    os.mkdir(log_path)
    for protein in proteins:
        os.mkdir(Path(output_path + '/' + protein))
        os.mkdir(Path(log_path + '/' + protein))
        trj = Path(trajectory_path + '/' + protein.split("_")[1] + '/a_input')
        bd = Path(boundary_path + '/' + protein.split("_")[1] + '/Center.bd')
        twn = Path(path + '/TWN/' + protein)
        out = Path(output_path + '/' + protein)
        log = Path(log_path + '/' + protein)
        run = f"-trj {trj.as_posix()} -bd {bd.as_posix()} -twn {twn.as_posix()} -o {out.as_posix()} -l {log.as_posix()}"
        subprocess.run(args=[sys.executable, 'TWN_Pattern.py'] + run.split(' '))

    # Region identification
    run = f"-d {Path(path).as_posix()}"
    subprocess.run(args=[sys.executable, 'TWN_Region.py'] + run.split(' '))