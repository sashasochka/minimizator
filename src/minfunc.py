#!/usr/bin/python3
import sys

verbose_messages = False
raw_result = False
helpmsg = (
    '''Usage: minfunc [options] [file]
     If file not specified read from stdin

     Command line options:
      -h [--help] for this help;
      -v [--verbose] to see additional minimization information.
      -r [--raw] to receive raw results for parsing. Not conformant with --verbose
      -n [--nohelp] to disable help at startup
    ''')

greetings = (
    '''This is Function Minimizator from Sochka Oleksandr!
    (Run with --help option to see additional info)
    Enter a system of your functions, each function on the separate line
    Enter function as a list of constituent numbers
    Press enter to finish''')


def main():
    """ Parses command line and interacts with user """
    global verbose_messages, raw_result
    nogreetings = False
    fin = sys.stdin
    for arg in sys.argv[1:]:
        if arg in ['-h', '--help']:
            print(helpmsg, end='')
            sys.exit(0)
        elif arg in ['-n', '--nohelp']:
            nogreetings = True;
        elif arg in ['-v', '--verbose']:
            verbose_messages = True
        elif arg in ['-r', '--raw']:
            raw_result = True
            nogreetings = True
        elif arg[0] == '-':
            sys.stderr.write('Option "{}" not recognized\n'.format(arg))
        else:
            try:
                fin = open(arg)
                nogreetings = True
            except IOError:
                sys.stderr.write('{} (File or folder not found!\n'.format(arg))
                sys.exit(2)

    if not nogreetings:
        print(greetings)

    system = []
    function_number = 1
    while True:
        if not nogreetings:
            print('Function {}: '.format(function_number), end='')
            sys.stdout.flush()
        line = fin.readline()
        if line is None or line.strip() == '':
            break
        if ';' in line: # has undefined
            posValsStr = line[: line.index(';')].strip()
            undefValsStr = line[line.index(';') + 1:].strip()
            undefVals = [int(i) for i in undefValsStr.split()]
        else:
            posValsStr = line
            undefVals = []

        posVals = [int(i) for i in posValsStr.split()]
        system.append((list(set(posVals + undefVals)), posVals))
        function_number += 1

    result = minimized(system)
    if not raw_result:
        print('Resulting minimal DNF system: ')
        for i, function in enumerate(result):
            print('Function {}: {}'.format(i, ' v '.join(function)))
    else:
        for function in result:
            print(' '.join(function))


def minimized(system):
    """ Minimizes the system. Accepts two-dimensional integer list """
    for function in system:
        for t in function:
            t = list(sorted(t))

    n = len(system) # the number of functions in the system
    k = get_number_of_args(system) # the number of arguments in the system
    verbose('\nThe number of arguments is {}...'.format(k))

    ''' calculate Full DNF'''
    phi_fdnf = [] # dnf of conjunction of the functions in the system
    values = list(value
                  for function in system
                  for t in function
                  for value in t)
    values = list(set(values))
    for v in values: # iterate over every possible value
        labels = []
        for i, func in enumerate(system):
            if binary_search(func[0], v) != -1:
                labels.append(i)
        if labels:
            phi_fdnf.append((bin(v, k), tuple(labels)))
    verbose('\n\nFull DNF:')
    verbose('\n'.join(map(str, phi_fdnf)))

    ''' calculate Reduced DNF '''
    phi_rdnf = []
    cur = phi_fdnf
    for i in range(k):
        cur_left = [True] * len(cur)
        next = []
        for a in range(len(cur)):
            for b in range(a + 1, len(cur)):
                common = tuple(set(cur[a][1]) & set(cur[b][1]))
                if not common:
                    continue
                diff = bit_diff(cur[a][0], cur[b][0])
                if diff == 0:
                    cur_left[a] = False
                elif diff == 1:
                    next.append((
                        combine(cur[a][0], cur[b][0]),
                        common
                    ))
                    if common == cur[a][1]:
                        cur_left[a] = False
                    if common == cur[b][1]:
                        cur_left[b] = False
        phi_rdnf.extend([cur[i] for i in range(len(cur)) if cur_left[i]])
        cur = unique(next)
    phi_rdnf.extend(cur)
    verbose('\n\nReduced DNF:')
    verbose('\n'.join(map(str, phi_rdnf)))

    ''' calculate Deadlock DNF using coverage table '''
    # <superhack>
    phi_fdnf_without_undef = [] # dnf of conjunction of the functions in the system
    values = list([value
                  for function in system
                  for value in function[1]])
    values = list(set(values))
    for v in values: # iterate over every possible value
        labels = []
        for i, func in enumerate(system):
            if binary_search(func[1], v) != -1:
                labels.append(i)
        if labels:
            phi_fdnf_without_undef.append((bin(v, k), tuple(labels)))
    # </superhack>

    columns = []
    for constituent, labels in phi_fdnf_without_undef:
        columns.extend([(constituent, label) for label in labels])
    rows = phi_rdnf

    # assert covers('XXX', '111')
    # assert covers('11X', '111')
    table = [[' '] * len(columns) for i in range(len(rows))]
    for i, (implicant, labels) in enumerate(rows):
        for j, (constituent, label) in enumerate(columns):
            if label in labels and covers(implicant, constituent):
                table[i][j] = 'Y'

    kernel_indexes = []
    for col in range(len(table[0])):
        index = -1
        for row in range(len(table)):
            if table[row][col] == 'Y':
                if index == -1:
                    index = row
                else:
                    break
        else:
            assert index != -1
            kernel_indexes.append(index)
    kernel_indexes = sorted(kernel_indexes)
    kernel = [rows[index] for index in kernel_indexes]
    verbose('\n\nKernel:')
    verbose('\n'.join(map(str, kernel)))


    non_kernel_cols = set(range(len(columns)))
    for i in kernel_indexes:
        for c in range(len(columns)):
            if table[i][c] == 'Y':
                if c in non_kernel_cols:
                    non_kernel_cols.remove(c)
    non_kernel_cols = list(non_kernel_cols)

    nonkernel = []
    l = 0
    for val in kernel_indexes:
        nonkernel.extend(range(l, val))
        l = val + 1
    nonkernel.extend(range(l, len(rows)))


    ''' calculate minimal DNF '''
    phi_mdnf = phi_rdnf
    # indexes = range(len(phi_rdnf))
    for mask in range((1 << len(nonkernel)) - 1, -1, -1):
        additional_indexes = [nonkernel[i]
                              for i, val in enumerate(bin(mask, len(nonkernel))) if val == '1']
        all_indexes = set(kernel_indexes) | set(additional_indexes)
        implicants = kernel + [rows[i] for i in additional_indexes]
        if complexity(implicants) >= complexity(phi_mdnf):
            continue

        for col in non_kernel_cols:
            local_ok = False
            for row in all_indexes:
                if table[row][col] == 'Y':
                    local_ok = True
                    break
            if not local_ok:
                break
        else:
            phi_mdnf = implicants
            # indexes = all_indexes

    verbose('\n\nMinimal DNF of phi:')
    verbose('\n'.join(map(str, phi_mdnf)))
    verbose()

    result = [set() for i in range(n)]
    for val, labels in phi_mdnf:
        for label in labels:
            result[label].add(val)

    for res in result:
        rmls = []
        for val1 in res:
            for val2 in res:
                if val1 != val2 and covers(val1, val2):
                    rmls.append(val2)
        for v in rmls:
            if v in res:
                res.remove(v)

    return result

def get_number_of_args(system):
    """ returns the maximal number of arguments from function of the system"""
    k = max([0] + [func[0][-1].bit_length() for func in system])
    return k

def bin(val, k):
    ''' converts val to binary with leading zeros such that it's length is equal to k '''
    return '{:0{}b}'.format(val, k)

def binary_search(a, val, lo=0, hi=None):
    ''' Performs binary search on sorted array and returns index if found or -1 else '''
    if hi is None:
        hi = len(a)
    while lo < hi:
        mid = (lo + hi) // 2
        midval = a[mid]
        if midval < val:
            lo = mid + 1
        elif midval > val:
            hi = mid
        else:
            return mid
    return -1

def unique(sequence):
    return list(set(sequence))

def bit_diff(a, b):
    """ returns the number of different bits on same positions in a and b """
    assert (len(a) == len(b))
    return sum(int(v1 != v2) for v1, v2 in zip(a, b))

def combine(a, b):
    """ aggregates a and b string into 1 string with different letters changed to X """
    assert (len(a) == len(b))
    return ''.join([v1 if v1 == v2 else 'X' for v1, v2 in zip(a, b)])

def covers(implicant, constituent):
    assert len(constituent) == len(implicant)
    return all(b in [a, 'X'] for a, b in zip(constituent, implicant))

def verbose(msg=''):
    global verbose_messages
    if verbose_messages:
        print(msg)

def complexity(implicants):
    result = 0
    for impl in implicants:
        result += len(impl[0]) - impl[0].count('X')
    return result - len(implicants)

if __name__ == '__main__':
    main()