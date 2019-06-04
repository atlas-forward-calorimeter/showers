"""Define command line run behavior.

Written by Anson Kost with the help of Professor John Rutherfoord, May 2019.

"""

import sys

if __name__ == '__main__':
    from core import analysis1

    if (len(sys.argv) == 1):
        # Testing only.
        print('Running from command line without arguments. This is for '
              'testing only.')
        analysis1('data')

    else:
        # Analyze given folder.
        assert(len(sys.argv) == 2)
        analysis1(sys.argv[1])
