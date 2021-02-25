#!/usr/bin/python3
"""
Generate rollsum chunking modelling graphs.

Usage: %(cmd)s [dir]
"""
from chunker import *
import pickle
import sys
import matplotlib
matplotlib.use('svg')
import matplotlib.pyplot as plt
from matplotlib.ticker import *

# Dimensions for graphs;
#
# chunker (chunker, Weibull1, Weibull2, FastCDCChunker)
# min len (0, 1/4, 1/2, 3/4)
# max_len (16, 8, 4, 2)
# avg_len (1, 2, 4, 8, 16, 32, 64)
#
# avg_size, max_size, min_size, dup_pct
#
# per chunker,
#   pct for min_len per avg_len, max_len=16
#   pct for max_len per avg_len, min_len=0
#   max_size for min_len per avg_len max_len=16.
#
# pct for avg_size by chunker, min_len=0, max_len=16
# pct for avg_size by chunker, min_len=1/2, max_len=2
# pct for avg_size by chunker, min_len=opt, max_len=opt(16?)


def FileName(dir,stat,alg,avg,min,max):
  """get svg graph file name from dimensions.

  For stat use 'perf|size(min|avg|max)' used for the y axis.
  Use 't' for the dimension used as the timeseries label.
  Use 'x' for dimension used as x axis.
  Use 'o' for the optimal value for the algorithm.
  Use 's' for the standard value for the algorithm.
  """
  #print("%r %r %r %r" % (alg, avg, min, max))
  if min not in tuple('txos'):
    min = '%.1f' % min
  if max not in tuple('txos'):
    max = '%.1f' % max
  return '%s/%s-%s-%s-%s-%s.svg' % (dir, stat, alg, avg, min, max)

def GetFileData(dir, alg):
  return pickle.load(open('%s/%s.dat' % (dir,alg), 'rb'))

def saveplt(filename, title, xlabel, ylabel, xticks=None, xlabels=None, yticks=None, ylabels=None):
  plt.title(title)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  if xticks:
    xlabels = xlabels or [str(n) for n in xticks]
    plt.xticks(xticks, xlabels)
  if yticks:
    ylabels = ylabels or [str(n) for n in yticks]
    plt.yticks(yticks, ylabels)
  #ax = plt.gca()
  #ax.ticklabel_format(axis='y', style='plain', useOffset=False)
  plt.grid()
  plt.legend(bbox_to_anchor=(1,1), loc='upper left')
  plt.savefig(filename, bbox_inches='tight')
  plt.cla()


def SizeMaxVsMinLimit(dir, data, avg, max):
  """Plot how size max varys with min limit."""
  d = data
  xs = mins
  for alg in algs:
    ys = [d[alg][avg][min][max][1].max/(avg * bsize) for min in mins]
    plt.plot(xs, ys, label="alg=%s" % alg)
  ax = plt.gca()
  ax.set_ylim(bottom=0)
  ax.set_xlim(left=0, right=xs[-1])
  #plt.xscale('log')
  #plt.yscale('log')
  saveplt(FileName(dir, 'sizemax','t', avg, 'x', max),
          'size max vs min limit for avg=%s,max=%.1fx' % (avg, max),
          'min limit fraction of avg', 'size max multiple of avg') #, xticks=xs)

def SizeDevVsMinLimit(dir, data, avg, max):
  """Plot how size dev varys with min limit."""
  d = data
  xs = mins
  for alg in algs:
    ys = [d[alg][avg][min][max][1].dev/(avg * bsize) for min in mins]
    plt.plot(xs, ys, label="alg=%s" % alg)
  ax = plt.gca()
  ax.set_ylim(bottom=0)
  ax.set_xlim(left=0, right=xs[-1])
  #plt.xscale('log')
  #plt.yscale('log')
  saveplt(FileName(dir, 'sizedev', 't', avg, 'x', max),
          'size dev vs min limit for avg=%s,max=%.1fx' % (avg, max),
          'min limit fraction of avg', 'size dev multiple of avg') #, xticks=xs)

def SizeAvgVsAvgTarget(dir, data, alg, max):
  """Plot how size avg varys with avg target."""
  d = data
  xs = avgs
  for min in mins:
    ys = [d[alg][avg][min][max][1].avg/(avg * bsize) for avg in avgs]
    plt.plot(xs, ys, label="min=%.1fx" % min)
  ax = plt.gca()
  #ax.set_ylim(bottom=0)
  #ax.set_xlim(left=0, right=xs[-1])
  plt.xscale('log')
  #plt.yscale('log')
  saveplt(FileName(dir, 'sizeavg', alg, 'x', 't', max),
          'size avg vs avg target for alg=%s,max=%.1fx' % (alg,max),
          'avg target', 'size avg multiple of avg', xticks=xs)

def PerfVsMinLimitByAvg(dir, data, alg, max):
  """Plot how deduplication peformance varys with min limit."""
  d = data[alg]
  xs = mins
  for avg in avgs:
    ys = [d[avg][min][max][0]*100.0 for min in mins]
    plt.plot(xs, ys, label="avg=%s" % avg)
  ax = plt.gca()
  ax.set_ylim(bottom=0)
  ax.set_xlim(left=0, right=xs[-1])
  #plt.xscale('log')
  #plt.yscale('log')
  saveplt(FileName(dir, 'perf', alg, 't', 'x', max),
          'Deduplication vs min limit for alg=%s,max=%.1fx' % (alg,max),
          'min limit fraction of avg', 'found %') #, xticks=xs)

def PerfVsMaxLimitByAvg(dir, data, alg, min):
  """Plot how deduplication peformance varys with max limit."""
  d = data[alg]
  xs = maxs
  for avg in avgs:
    ys = [d[avg][min][max][0]*100.0 for max in maxs]
    plt.plot(xs, ys, label="avg=%s" % avg)
  ax = plt.gca()
  ax.set_ylim(bottom=0)
  #ax.set_xlim(left=0, right=xs[-1])
  plt.xscale('log')
  #plt.yscale('log')
  saveplt(FileName(dir, 'perf', alg, 't', min, 'x'),
          'Deduplication vs max limit for alg=%s,min=%.1fx' % (alg,min),
          'max limit multiple of avg', 'found %', xticks=xs)

def PerfVsMinLimitByAlg(dir, data, avg, max):
  """Plot how deduplication peformance varys with min limit."""
  d = data
  xs = mins
  for alg in algs:
    ys = [d[alg][avg][min][max][0]*100.0 for min in mins]
    plt.plot(xs, ys, label="alg=%s" % alg)
  ax = plt.gca()
  ax.set_ylim(bottom=0)
  ax.set_xlim(left=0, right=xs[-1])
  #plt.xscale('log')
  #plt.yscale('log')
  saveplt(FileName(dir, 'perf', 't', avg, 'x', max),
          'Deduplication vs min limit for avg=%s,max=%.1fx' % (avg,max),
          'min limit fraction of avg', 'found %') #, xticks=xs)

def PerfVsAvgSize(dir, data, min, max):
  """Plot how deduplication peformance varys with avg size."""
  d = data
  for alg in algs:
    xs = [d[alg][avg][min][max][1].avg/bsize for avg in avgs]
    ys = [d[alg][avg][min][max][0]*100.0 for avg in avgs]
    plt.plot(xs, ys, label="alg=%s" % alg)
  ax = plt.gca()
  ax.set_ylim(bottom=0)
  #ax.set_xlim(left=0, right=xs[-1])
  plt.xscale('log')
  #plt.yscale('log')
  saveplt(FileName(dir, 'perf', 't', 'x', min, max),
          'Deduplication vs avg target for min=%.1fx,max=%.1fx' % (min, max),
          'avg size', 'found %', xticks=avgs)

dir = sys.argv[1] if len(sys.argv) > 1 else '.'
if dir in ("-?", "-h", "--help", None):
  print(__doc__ % dict(cmd=os.path.basename(sys.argv[0])))
  sys.exit(1)

# This builds a results dict structured as;
# results[alg][avg][min][max] = (perf, blkstats, dupstats)
#algs = sorted(chunkers)
algs = 'chunker weibull1 weibull2 weibullt1 weibullt2 nc1 nc2 nc3 rc4'.split()
data = {}
for alg in algs:
  tsize, bsize, data[alg] = GetFileData(dir, alg)

alg0 = data[algs[0]]
avgs = sorted(alg0)
avg0 = alg0[avgs[0]]
mins = sorted(avg0)
min0 = avg0[mins[0]]
maxs = sorted(min0)

for alg in algs:
  PerfVsMinLimitByAvg(dir, data, alg, maxs[0])
  PerfVsMinLimitByAvg(dir, data, alg, maxs[-1])
  PerfVsMaxLimitByAvg(dir, data, alg, 0.0)
  PerfVsMaxLimitByAvg(dir, data, alg, 0.5)
  SizeAvgVsAvgTarget(dir, data, alg, maxs[0])
  SizeAvgVsAvgTarget(dir, data, alg, maxs[-1])
#SizeMaxVsMinLimit(dir, data, avgs[0], maxs[-1])
SizeDevVsMinLimit(dir, data, avgs[0], maxs[-1])
for avg in avgs:
  PerfVsMinLimitByAlg(dir, data, avg, maxs[0])
  PerfVsMinLimitByAlg(dir, data, avg, maxs[-1])
for min in mins:
  PerfVsAvgSize(dir, data, min, maxs[0])
  PerfVsAvgSize(dir, data, min, maxs[-1])
PerfVsAvgSize(dir, data, 0.5, 2.0)
PerfVsAvgSize(dir, data, 0.5, 4.0)
