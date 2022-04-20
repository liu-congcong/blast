#!/usr/bin/env python3
import os
import sys
from argparse import ArgumentParser, RawTextHelpFormatter
from datetime import datetime
from math import ceil
from multiprocessing import Process
from shutil import which
from subprocess import DEVNULL, run
from tempfile import mkstemp


def __init__(parameters):
    parser = ArgumentParser(
        formatter_class = RawTextHelpFormatter,
        description = 'Run blast with multiple threads.',
        epilog = 'Author: Liucongcong.',
    )

    parser.add_argument(
        '-q', '--query', type = str, required = True, metavar = 'fasta'
    )

    parser.add_argument(
        '-t', '--target', type = str, required = True, metavar = 'database|fasta'
    )

    parser.add_argument(
        '-f', '--function', type = str, required = True, choices=('blastn', 'blastp', 'blastx', 'tblastn', 'tblastx'),
        metavar = 'blastn|blastp|blastx|tblastn|tblastx'
    )

    parser.add_argument(
        '-o', '--output', type = str, required = True, metavar = 'output',
    )

    parser.add_argument(
        '--makeblastdb', default = which('makeblastdb', mode = os.X_OK), type = str, required = False, metavar = 'makeblastdb',
        help = 'Default: environmental makeblastdb.'
    )

    parser.add_argument(
        '--blastn', default = which('blastn', mode = os.X_OK), type = str, required = False, metavar = 'blastn',
        help = 'Default: environmental blastn.'
    )

    parser.add_argument(
        '--blastp', default = which('blastp', mode = os.X_OK), type = str, required = False, metavar = 'blastp',
        help = 'Default: environmental blastp.'
    )

    parser.add_argument(
        '--blastx', default = which('blastx', mode = os.X_OK), type = str, required = False, metavar = 'blastx',
        help = 'Default: environmental blastx.'
    )

    parser.add_argument(
        '--tblastn', default = which('tblastn', mode = os.X_OK), type = str, required = False, metavar = 'tblastn',
        help = 'Default: environmental tblastn.'
    )

    parser.add_argument(
        '--tblastx', default = which('tblastx', mode = os.X_OK), type = str, required = False, metavar = 'tblastx',
        help = 'Default: environmental tblastx.'
    )

    parser.add_argument(
        '--threads', type = int, default = os.cpu_count(), required = False, metavar = 'threads',
        help = 'Default: {0}.'.format(os.cpu_count())
    )

    parser.add_argument(
        '--strand', type = str, default = 'both', required = False, choices = ('both', 'minus', 'plus'), metavar = 'both|minus|plus',
        help = 'Default: both.'
    )

    parser.add_argument(
        '--codon_table', type = int, default = 1, required = False, metavar = '1-6, 9-16, 21-25',
        choices=(1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25),
        help = '1: The Standard Code\n2: The Vertebrate Mitochondrial Code\n3: The Yeast Mitochondrial Code\n4: The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code\n5: The Invertebrate Mitochondrial Code\n6: The Ciliate, Dasycladacean and Hexamita Nuclear Code\n9: The Echinoderm and Flatworm Mitochondrial Code\n10: The Euplotid Nuclear Code\n11: The Bacterial, Archaeal and Plant Plastid Code\n12: The Alternative Yeast Nuclear Code\n13: The Ascidian Mitochondrial Code\n14: The Alternative Flatworm Mitochondrial Code\n15: ??? Code\n16: Chlorophycean Mitochondrial Code\n21: Trematode Mitochondrial Code\n22: Scenedesmus obliquus Mitochondrial Code\n23: Thraustochytrium Mitochondrial Code\n24: Pterobranchia Mitochondrial Code\n25: Candidate Division SR1 and Gracilibacteria Code\nDefault: 1.'
    )
    return parser.parse_args(parameters)


def make_file():
    file_descriptor, file_name = mkstemp(dir = os.getcwd())
    os.close(file_descriptor)
    return file_name


def run_makeblastdb(makeblastdb, dbtype, input_fasta):
    output_blastdb = make_file()
    run_process = run(
        [makeblastdb, '-in', input_fasta, '-dbtype', dbtype, '-hash_index', '-out', output_blastdb],
        stdout = DEVNULL, stderr = DEVNULL
    )
    assert not run_process.returncode, 'An error has occured while running makeblastdb.'
    return output_blastdb


def split_fasta(input_file, n):
    position_list = list()
    open_input = open(input_file, 'rb')
    while True:
        line = open_input.readline()
        if line.startswith(b'>'):
            position_list.append(open_input.tell() - len(line))
        if not line:
            position_list.append(open_input.tell())
            break
    step = ceil((len(position_list) - 1) / n)
    position_list_index = 0
    while True:
        positions = position_list[position_list_index : position_list_index + step + 1]
        if len(positions) <= 1:
            break
        output_file = make_file()
        open_input.seek(positions[0], 0)
        open_output = open(output_file, 'wb')
        open_output.write(open_input.read(positions[-1] - positions[0]))
        open_output.close()
        yield (output_file)
        position_list_index += step
    open_input.close()
    return None


def run_blast_thread(command, input_file, output_file):
    run_process = run(
        command + ['-query', input_file, '-out', output_file],
        stdout = DEVNULL, stderr = DEVNULL
    )
    assert not run_process.returncode, 'An error has occured while running blast.'
    os.remove(input_file)
    return None


def remove_blastdb(blastdb_prefix):
    open_directory = os.scandir(path = os.getcwd())
    for entry in open_directory:
        if entry.is_file() and entry.path.startswith(blastdb_prefix):
            os.remove(entry.path)
    open_directory.close()
    return None


def combine(input_files, output_file):
    open_output = open(output_file, 'wb')
    header = '\t'.join(['qseqid', 'qstart', 'qend', 'qlen', 'sseqid', 'sstart', 'send', 'slen', 'pident', 'score'])
    open_output.write(header.encode() + b'\n')
    for input_file in input_files:
        open_input = open(input_file, 'rb')
        while True:
            if not open_output.write(open_input.read(1024 ** 3)):
                break
        open_input.close()
        os.remove(input_file)
    open_output.close()
    return None


if '__main__' == __name__:
    parameters = __init__(sys.argv[1 : ])

    assert parameters.makeblastdb and os.access(parameters.makeblastdb, os.X_OK), '"--makeblastdb" should be specified.'
    assert parameters.blastn and os.access(parameters.blastn, os.X_OK), '"--blastn" should be specified.'
    assert parameters.blastp and os.access(parameters.blastp, os.X_OK), '"--blastp" should be specified.'
    assert parameters.blastx and os.access(parameters.blastx, os.X_OK), '"--blastx" should be specified.'
    assert parameters.tblastn and os.access(parameters.tblastn, os.X_OK), '"--tblastn" should be specified.'
    assert parameters.tblastx and os.access(parameters.tblastx, os.X_OK), '"--tblastx" should be specified.'

    if os.access(parameters.target, os.R_OK): # fasta file #
        makeblastdb_marker = True
        dbtype = r'nucl' if parameters.function in (r'blastn', r'tblastn', r'tblastx') else r'prot'
        print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '->', 'Making database for blast.', flush = True)
        parameters.target = run_makeblastdb(parameters.makeblastdb, dbtype, parameters.target)
        print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '->', 'Done.', flush = True)
    else:
        makeblastdb_marker = False

    command = [
        parameters.function, '-db', parameters.target, '-num_threads', '1',
        '-outfmt', '6 qseqid qstart qend qlen sseqid sstart send slen pident score',
    ]

    # construct blast command #
    if parameters.function == 'blastn':
        command.extend(['-strand', parameters.strand])
    elif parameters.function == 'blastx':
        command.extend(['-strand', parameters.strand, '-query_gencode', str(parameters.codon_table)])
    elif parameters.function == 'tblastx':
        command.extend(['-strand', parameters.strand, '-db_gencode', str(parameters.codon_table), '-query_gencode', str(parameters.codon_table)])
    elif parameters['function'] == 'tblastn':
        command.extend(['-db_gencode', str(parameters.codon_table)])

    # run blast #
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '->', 'Running blast.', flush = True)
    process_list = list()
    output_list = list()
    for query_file in split_fasta(parameters.query, parameters.threads):
        output_list.append(make_file())
        process_list.append(
            Process(
                target = run_blast_thread,
                args = (command, query_file, output_list[-1]),
            )
        )
    for process in process_list:
        process.start()
    for process in process_list:
        process.join()
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '->', 'Done.', flush = True)

    if makeblastdb_marker:
        remove_blastdb(parameters.target)

    combine(output_list, parameters.output)

    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '->', 'Finished.', flush = True)
