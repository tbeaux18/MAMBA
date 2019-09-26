#!/usr/bin/env python3
"""
@author: Timothy Baker

fastqc.py

- will run fastqc at 24 threads
"""

import pathlib
import argparse
import shlex
import subprocess
import logging
import sys
import re
import collections


logger = logging.getLogger(__name__)


def arg_parser():
    """ Argument input from command line """

    parser = argparse.ArgumentParser(
        description='Script takes a directory of fastq.gz files and runs FASTQC'
    )
    parser.add_argument('-d', '--fastq-dir', type=str, help='Enter path/to/dir_of_fastqs')

    return parser.parse_args()


def gather_fastq(fastq_directory):
    """ parses the fastq directory, parses names and aligns read 1 and read 2 """

    # regex specific for this project directory
    sample_regex = re.compile(r'(.*?)_L001')

    fastq_dict = collections.defaultdict(list)

    for path in fastq_directory.iterdir():

        sample_name = sample_regex.match(path.stem)

        if sample_name:
            fastq_dict[sample_name.group(1)].append(str(path))
        else:
            print(path, "is probably a directory.")

    return fastq_dict


def run_fastqc(fastqc_log, **kwargs):
    """ runs fastqc with the provided args
        params:
            kwargs : dict specified in build_fastqc_args
            --extract : qc log files will be unzipped
            --format : fastq default
        returns:
            None
    """

    fastqc_command = """fastqc
                        --extract
                        --threads 20
                        --format fastq
                        --outdir {output_dir}
                        {fastq_read1}
                        {fastq_read2}""".format(**kwargs)

    fastqc_formatted_args = shlex.split(fastqc_command)

    with open(str(fastqc_log), 'ab+') as fastqc_stdout:
        fastqc_process = subprocess.Popen(
            fastqc_formatted_args,
            bufsize=20,
            stdout=fastqc_stdout,
            stderr=fastqc_stdout
        )
        fastqc_process.communicate()



def main():
    """ runs main """

    args = arg_parser()

    # only pass in directory of demultiplexed fastq
    if args.fastq_dir:
        fastq_dir = pathlib.Path(args.fastq_dir)

    fastq_main_dict = gather_fastq(fastq_dir)

    script_path = pathlib.Path(__file__)

    script_parent = script_path.resolve().parent

    log_path = script_parent.joinpath('logs')

    log_path.mkdir(exist_ok=True)

    log_file_path = log_path.joinpath('fastqc_log.txt')

    log_file_path.touch(exist_ok=True)

    for sam_name, fastq_array in fastq_main_dict.items():

        # creating a sub directory for fastqc output
        fastq_output = fastq_dir.joinpath(sam_name)

        # creating output directory before passing to fastqc
        fastq_output.mkdir(exist_ok=True)

        run_fastqc(
            log_file_path,
            output_dir=str(fastq_output),
            fastq_read1=fastq_array[0],
            fastq_read2=fastq_array[1]
        )



if __name__ == '__main__':
    main()
