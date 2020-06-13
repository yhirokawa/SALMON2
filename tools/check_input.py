#!/usr/bin/env python
#
#   Copyright 2020 ARTED developers
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#

import math
import re
import sys

def list_div(l1,l2):
  n = min(len(l1),len(l2))
  p = []
  for i in range(0,n):
    p.append(l1[i] / l2[i])
  return p

def list_prod(nlist):
  n = 1
  for i in nlist:
    n = n * i
  return n

def fail_exit():
  print 'ERROR.'
  sys.exit(1)

# SALMON input list
def to_boolean(s):
  if s.lower() == '\'y\'':
    return True
  elif s.lower() == '\'n\'':
    return False
  else:
    return False

def gen_inputlist_map(filename):
  inputmap = {}
  inputmap['nproc_k']            = 1
  inputmap['nproc_ob']           = 1
  inputmap['nproc_rgrid']        = [1,1,1]
  inputmap['yn_ffte']            = False
  inputmap['process_allocation'] = 'grid_sequential'
  with open(filename) as f:
    number_pattern      = re.compile(r'\d+')
    float_pattern       = re.compile(r'[+-]?(?:\d+\.?\d*|\.\d+)(?:[dD][+-]?\d+)?')
    yn_pattern          = re.compile(r'\'[yYnN]\'')
    string_pattern      = re.compile(r'\'\w*\'')
    theory_pattern      = re.compile(r'\s*theory\s*=\s*')
    nproc_k_pattern     = re.compile(r'\s*nproc_k\s*=\s*\d+')
    nproc_ob_pattern    = re.compile(r'\s*nproc_ob\s*=\s*\d+')
    nproc_rgrid_pattern = re.compile(r'\s*nproc_rgrid\s*=\s*\d+')
    al_pattern          = re.compile(r'\s*al\s*=\s*\d+')
    dl_pattern          = re.compile(r'\s*dl\s*=\s*\d+')
    num_rgrid_pattern   = re.compile(r'\s*num_rgrid\s*=\s*\d+')
    yn_ffte_pattern     = re.compile(r'\s*yn_ffte\s*=\s*')
    proc_alloc_pattern  = re.compile(r'\s*process_allocation\s*=\s*')
    for s in f:
      if theory_pattern.match(s):
        inputmap['theory'] = string_pattern.search(s).group().strip('\'').lower()
      elif nproc_k_pattern.match(s):
        inputmap['nproc_k'] = int(number_pattern.search(s).group())
      elif nproc_ob_pattern.match(s):
        inputmap['nproc_ob'] = int(number_pattern.search(s).group())
      elif nproc_rgrid_pattern.match(s):
        dims = []
        for x in number_pattern.findall(s):
          dims.append(int(x))
        inputmap['nproc_rgrid'] = dims
      elif al_pattern.match(s):
        dims = []
        for x in float_pattern.findall(s):
          dims.append(float(x.replace('d','e')))
        inputmap['al'] = dims
      elif dl_pattern.match(s):
        dims = []
        for x in float_pattern.findall(s):
          dims.append(float(x.replace('d','e')))
        inputmap['dl'] = dims
      elif num_rgrid_pattern.match(s):
        dims = []
        for x in number_pattern.findall(s):
          dims.append(int(x))
        inputmap['num_rgrid'] = dims
      elif yn_ffte_pattern.match(s):
        inputmap['yn_ffte'] = to_boolean(yn_pattern.search(s).group())
      elif proc_alloc_pattern.match(s):
        inputmap['process_allocation'] = string_pattern.search(s).group().strip('\'').lower()
  return inputmap

def print_inputlists(inputmap):
  print '&parallel'
  print '  nproc_k = {}'.format(inputmap['nproc_k'])
  print '  nproc_ob = {}'.format(inputmap['nproc_ob'])
  print '  nproc_rgrid = {},{},{}'.format(*inputmap['nproc_rgrid'])
  print '  process_allocation = \'{}\''.format(inputmap['process_allocation'])
  print '/'
  print '&rgrid'
  print '  num_rgrid = {},{},{}'.format(*inputmap['num_rgrid'])
  print '/'


# FFTE
def check_ffte_condition_gridsize(num_rgrid, nprocs_rgrid):
  y1 = num_rgrid[0] % nprocs_rgrid[1]
  y2 = num_rgrid[1] % nprocs_rgrid[1]
  z1 = num_rgrid[1] % nprocs_rgrid[2]
  z2 = num_rgrid[2] % nprocs_rgrid[2]
  return y1 == 0 and y2 == 0 and z1 == 0 and z2 == 0

def check_ffte_condition_prime_factors(num_rgrid):
  for i in range(0,3):
    t = num_rgrid[i]
    for j in range(0,26):
      if t % 2 == 0:
        t = t / 2
    for j in range(0,17):
      if t % 3 == 0:
        t = t / 3
    for j in range(0,11):
      if t % 5 == 0:
        t = t / 5
    if t != 1:
      return False
  return True


def get_nproc_rgrid(nprocs, nprocs_per_node, num_rgrid, is_ffte):
  if num_rgrid[0] % nprocs_per_node != 0:
    print '[INFO] num_rgrid[0] % nprocs_per_node is not divided.'
  nzy = nprocs / nprocs_per_node
  nz  = int(math.sqrt(float(nzy)))
  nz  = nz + (nz % 2)
  for ii in range(0,100):
    ny = nzy / nz
    if (ny*nz == nzy):
      break
    nz = nz + 2
  if is_ffte:
    ny = ny + (ny % 2)
    for ii in range(ny,1,-1):
      if check_ffte_condition_gridsize(num_rgrid, [nprocs_per_node, ii, nz]):
        if nprocs % list_prod([nprocs_per_node, ii, nz]) == 0:
          ny = ii
          break
  return [nprocs_per_node, ny, nz]


def prime_factorize(n):
  a = []
  while n % 2 == 0:
    a.append(2)
    n //= 2
  f = 3
  while f * f <= n:
    if n % f == 0:
      a.append(f)
      n //= f
    else:
      f += 2
  if n != 1:
    a.append(n)
  return a

def determine_node_shape(nproc_k, nproc_ob, nproc_rgrid, nnodes, nprocs_per_node, process_order):
  if process_order == 'grid_sequential':
    PS    = nproc_ob*nproc_k
    shape = [nproc_rgrid[0], nproc_rgrid[1], nproc_rgrid[2]]
  else:
    PS    = nproc_rgrid[0]
    shape = [nproc_ob*nproc_k, nproc_rgrid[1], nproc_rgrid[2]]

  # determine # of sub-process grid
  # PS = NS[0] * NS[1] * NS[2]
  if PS == 1:
    NS = [1, 1, 1]
  else:
    factor = prime_factorize(PS)
    

  shape[0] = shape[0] * NS[0] / nprocs_per_node
  shape[1] = shape[1] * NS[1]
  shape[2] = shape[2] * NS[2]

  return shape


if __name__ == '__main__':
  if len(sys.argv) < 4:
    print '[Usage] ./{} <SALMON inputfile> <required # of node> <required # of procs/node>'.format(sys.argv[0])
    fail_exit()

  inputmap = gen_inputlist_map(sys.argv[1])

  # logical checking...
  if not 'num_rgrid' in inputmap.keys():
    dims = list_div(inputmap['al'], inputmap['dl'])
    inputmap['num_rgrid'] = [int(f) for f in dims]
    print '[INFO] num_rgrid constructed = {}'.format(inputmap['num_rgrid'])


  nnodes          = int(sys.argv[2])
  nprocs_per_node = int(sys.argv[3])
  nprocs_global   = nnodes * nprocs_per_node
  if list_prod(inputmap['nproc_rgrid']) * inputmap['nproc_ob'] * inputmap['nproc_k'] != nprocs_global:
    print '[INFO]'
    print 'product(nproc_rgrid) * nproc_ob * nproc_k /= # of MPI procs = {}'.format(nprocs_global)
    print 'calculate nproc_k,ob and rgrid'

    # find nproc_.*
    num_rgrid   = inputmap['num_rgrid']
    nproc_rgrid = get_nproc_rgrid(nprocs_global, nprocs_per_node, num_rgrid, inputmap['yn_ffte'])
    nproc_ob    = nprocs_global / list_prod(nproc_rgrid)
    nproc_k     = 1
    inputmap['nproc_k']     = nproc_k
    inputmap['nproc_ob']    = nproc_ob
    inputmap['nproc_rgrid'] = nproc_rgrid

  theory = inputmap['theory']
  if theory == 'dft' or theory == 'dft_md':
    # DFT prefers orbital_sequential
    inputmap['process_allocation'] = 'orbital_sequential'
  else:
    # TDDFT prefers grid_sequential
    inputmap['process_allocation'] = 'grid_sequential'


  # FFTE checking...
  if inputmap['yn_ffte']:
    if check_ffte_condition_gridsize(inputmap['num_rgrid'], inputmap['nproc_rgrid']):
      print '[FFTE]'
      print 'num_rgrid and nproc_rgrid are available to use FFTE.'
      if not check_ffte_condition_prime_factors(inputmap['num_rgrid']):
        print '       prime factors for number of grids must be combination of 2,3, or 5'
    else:
      print '[FFTE]'
      print 'num_rgrid and nproc_rgrid are unsuitable for using FFTE.'
      print 'please check for condition:'
      print '  mod(num_rgrid(1),nprocs_rgrid(2)) must be 0'
      print '  mod(num_rgrid(2),nprocs_rgrid(2)) must be 0'
      print '  mod(num_rgrid(2),nprocs_rgrid(3)) must be 0'
      print '  mod(num_rgrid(3),nprocs_rgrid(3)) must be 0'
      fail_exit()


  print ''
  print '# =============================================== #'
  print 'Probably suitable parameters for large scale system.'
  print 'please replace the following inputlists.'
  print ''
  print_inputlists(inputmap)

  print ''
  print '# =============================================== #'
  print 'If you use Tofu-D network or other 3-D network.'
  print 'probably, the following node shape is suitable to run the application.'
  print ''
  print '  Requested network shape = (PX, PY, PZ, PW)'
  print '  Tofu-network shape      = (TX, TY, TZ)'
  print '  Process shape           = (TX*PPN, TY, TZ)'
  print ''
  print '  (PX, PY, PZ) = nproc_rgrid'
  print '  PW           = nproc_ob*nproc_k'
  print '  PPN          = # of process/node'
  if inputmap['process_allocation'] == 'grid_sequential':
    print '  PW  = PW1*PW2*PW3'
    print '  PW1 = TX*PPN / PX'
    print '  PW2 = TY     / PY'
    print '  PW3 = TZ     / PZ'
    print '  TX  = PX*PW1 / PPN'
    print '  TY  = PY*PW2'
    print '  TZ  = PZ*PW3'
  else:
    print '  PX  = PX1*PX2*PX3'
    print '  PX1 = TX*PPN / PW'
    print '  PX2 = TY     / PY'
    print '  PX3 = TZ     / PZ'
    print '  TX  = PW*PX1 / PPN'
    print '  TY  = PY*PX2'
    print '  TZ  = PZ*PX3'
  print ''
  print 'We determine the node shape to be as close to a regular hexahedron as possible.'
  shape = determine_node_shape(inputmap['nproc_k'], inputmap['nproc_ob'], inputmap['nproc_rgrid'] \
      , nnodes, nprocs_per_node, inputmap['process_allocation'])

  print '(TX, TY, TZ) = ({}, {}, {})'.format(shape[0], shape[1], shape[2])
