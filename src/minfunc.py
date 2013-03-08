#!/usr/bin/python3
import sys
import math
import itertools
import operator

verbose_messages = False
helpmsg = (
'''Usage: minfunc [options] [file]
 If file not specified read from stdin

 Command line options:
  -h [--help] for this help;
  -v [--verbose] to see additional minimization information
  -n [--nohelp] to disable help at startup
''')

greetings = (
'''This is Function Minimizator from Sochka Oleksandr!
(Run with --help option to see additional info)
Enter a system of your functions, each function on the separate line
Enter function as a list of constituent numbers
Press enter to finish''')

def main():
  ''' Parses command line and interacts with user '''
  global verbose_messages
  nogreetings = False
  fin = sys.stdin
  for arg in sys.argv[1:]:
    if arg in ['-h', '--help']:
      print(helpmsg, end = '')
      sys.exit(0)
    elif arg in ['-n', '--nohelp']:
      nogreetings = True;
    elif arg in ['-v', '--verbose']:
      verbose_messages = True
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
      print('Function {}: '.format(function_number), end = '')
      sys.stdout.flush()
    line = fin.readline()
    if line is None or line.strip() == '':
      break
    system.append(list(set(map(int, line.split()))))
    function_number += 1

  result = minimized(system)


def minimized(system):
  ''' Minimizes the system. Accepts two-dimensional integer list '''

  n = len(system) # the number of functions in the system
  k = get_number_of_args(system) # the number of arguments in the system
  print('The number of arguments is {}...'.format(k))

  ''' conjunct all functions in system to function phi to receive phi FDNF'''
  phi_fdnf = [] # dnf of conjunction of the functions in the system
  values = list({value for function in system for value in function })
  for i in values: # iterate over every possible value
    labels = []
    for j, func in enumerate(system):
      if binary_search(func, i) != -1:
        labels.append(j)
    if labels != []:
      phi_fdnf.append((bin(i,k), tuple(labels)))
  verbose('\n\nPhi DNF:')
  verbose('\n'.join(map(str, phi_fdnf)))

  ''' calculate SDNF '''
  phi_sdnf = []
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
    phi_sdnf.extend([cur[i] for i in range(len(cur)) if cur_left[i]])
    cur = unique(next)
  phi_sdnf.extend(cur)
  verbose('\n\nPhi SDNF:')
  verbose('\n'.join(map(str, phi_sdnf)))



def get_number_of_args(system):
  ''' returns the maximal number of arguments from function of the system'''
  k = 0
  for func in system:
    for val in func:
      if val != 0:
        k = max(k, int(math.log(val, 2.)) + 1)
  return k
  
def bin(val, k):
  ''' converts val to binary with leading zeros such that it length = k '''
  res = ['0'] * k
  i = k - 1
  while val != 0:
    res[i] = str(val%2)
    val //= 2 
    i -= 1
  return ''.join(res)

def binary_search(a, val, lo = 0, hi = None):
  ''' Performs binary search on sorted array and returns index if found or -1 else '''
  if hi is None:
      hi = len(a)
  while lo < hi:
      mid = (lo+hi)//2
      midval = a[mid]
      if midval < val:
          lo = mid+1
      elif midval > val: 
          hi = mid
      else:
          return mid
  return -1

def unique(sequence):
    return list(set(sequence))

def bit_diff(a, b):
  ''' returns the number of different bits on same positions in a and b '''
  assert(len(a) == len(b))
  result = 0
  for i in range(len(a)):
    if a[i] != b[i]:
      result += 1
  return result

def combine(a, b):
  ''' aggregates a and b string into 1 string with different letters changed to X '''
  assert(len(a) == len(b))
  return ''.join([v1 if v1 == v2 else 'X' for v1, v2 in zip(a, b)])

def verbose(msg):
  global verbose_messages
  if verbose_messages:
    print(msg)

if __name__ == '__main__':
  main()