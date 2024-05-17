"""!
@brief Run scripts/snp_analyser.py.
or any information refer to README.md or  scripts/snp_analyser.py.
"""
from scripts import snp_analyser as snp
import sys

if __name__ == '__main__':
    snp.main_using_getopts(sys.argv[1:])

