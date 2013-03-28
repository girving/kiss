#!/usr/bin/env python

from __future__ import division
from other.core import *
from kiss import *

def test_paths():
  # Known information about shortest paths in A_12
  exact = asarray([1,48,2016,80700,2891295,73385595,163078590,62495,60])
  assert sum(exact)==239500800
  exact_density = exact/sum(exact)
  print 'exact density = %s'%list(exact_density)
  # Evaluate various shortest paths
  counts = zeros(9,dtype=int) 
  random.seed(239500800)
  print 'generators = %s'%' '.join(map(cycles,generators()))
  for _ in xrange(2000):
    g = arange(12,dtype=int32)
    random.shuffle(g)
    g = from_permutation(g)
    if s12_parity(g)>0: 
      path = shortest_path(g)
      counts[len(path)] += 1
      if 0:
        print '%s = %d'%(cycles(g),len(path))
      assert reduce(s12_times,path,s12_identity)==g
  # Compare densities
  rough_density = counts/sum(counts)
  print 'rough density = %s'%list(rough_density)
  error = maxabs(rough_density-exact_density)
  print 'density error = %g'%error
  assert error<.01

if __name__=='__main__':
  test_paths()
