#!/usr/bin/env python

#import pylab as pl


import numarray as na
import numarray.mlab as mlab
import re

import pylab as pl
#import matplotlib as mplib
#import matplotlib.matlab as pl
#import matplotlib as mpl
#import matplotlib.pylab as pl

def file2na( path ):
  f = open(path)
  lines = f.readlines()
  lines = [ l for l in lines if re.match("^Trial:", l) ]
  rows = len(lines)
  if rows < 1:
    print "No rows dummy!"
    sys.exit(1)
  cols = len(re.findall('[0-9.e+-]+',lines[1]))
  a = na.array( shape = ( rows, cols ), type = na.Float64 )
  for i in xrange( rows ):
    line = re.findall('[0-9.e+-]+',lines[i])
    if len(line) != cols:
      nline = [ '0.0' for j in xrange(cols) ]
      for c in xrange(len(line)):
        nline[c] = line[c]
      line = nline
    row = [ eval(s) for s in line ]
    a[ i ] = row
  return a

if __name__ == "__main__":
  import sys
  import glob
  files = glob.glob(sys.argv[1])
  column_to_plot = 1 # begins with 1
  first_row_to_plot = 0 # to produce more detailed graphs... 15 for nac in sinus,
  i = 0
  b = na.zeros( shape = (1, len(files)), type = na.Float64 )
  for f in files:
    print f
    a = file2na( f )
    if i == 0:
      b = na.zeros( shape = (a.shape[0], len(files)), type = na.Float64 )
    b[:,i] += a[:,column_to_plot]
    i += 1
  last_row_to_plot = a.shape[0] # 25/30 for pg, 99 for nac
  print b.shape
  m = mlab.mean(b[first_row_to_plot:last_row_to_plot,:], 1)
  s = mlab.std(b[first_row_to_plot:last_row_to_plot,:], 1)
  min = mlab.min(b[first_row_to_plot:last_row_to_plot,:], 1)
  max = mlab.max(b[first_row_to_plot:last_row_to_plot,:], 1)
  if a[0,0] == 0.0:
    # Can't take log of 0. 
    a[0,0] = 1.0;
  
  lw = 2
  #pl.figure(figsize=(4, 3), dpi=300)
  pl.hold(True)

  pl.plot(a[first_row_to_plot:last_row_to_plot,0],m,'b-', label="Mean", linewidth=lw)
  pl.plot(a[first_row_to_plot:last_row_to_plot,0],m+s,'c--', label="Upper sdev.", linewidth=lw)
  pl.plot(a[first_row_to_plot:last_row_to_plot,0],m-s,'c--', label="Lower sdev.", linewidth=lw)
  pl.plot(a[first_row_to_plot:last_row_to_plot,0],min,'g:', label="Min", linewidth=lw)
  pl.plot(a[first_row_to_plot:last_row_to_plot,0],max,'r:', label="Max", linewidth=lw)

  pl.ylabel("R") # Reward / performance criteria
  pl.xlabel("Steps")

  #pl.axis([0,2000000,10,100])
  pl.legend(loc='lower right')

  pl.grid(True)
  pl.show()


