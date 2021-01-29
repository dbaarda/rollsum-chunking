#!/usr/bin/python3
"""
RollsumChunking modelling graphs.
"""
from chunker import *
import pickle
import matplotlib
matplotlib.use('svg')
import matplotlib.pyplot as plt
from matplotlib.ticker import *

# Dimensions for graphs;
#
# chunker (chunker, Weibull1, Weibull2, FastCDCChunker)
# min len (0, 1/4, 1/2, 3/4)
# max_len (16, 8, 4, 2)
# avg_len (1,2,4, 8, 16, 32, 64)
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


def FileName(stat,alg,min,avg,max):
  """get svg graph file name from dimensions.

  For stat use 'perf|size(min|avg|max)' used for the y axis.
  Use 'l' for the dimension used as the timeseries label.
  Use 'x' for dimension used as x axis.
  """
  return 'data/%s-%s-%s-%s-%s.svg' % (stat, alg, min, avg, max)

def GetFileData(alg):
  return pickle.load(open('data/%s.dat' % alg, 'rb'))

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


def SizeMaxVsMinLimit(data, avg, max):
  """Plot how size max varys with min limit."""
  d = data
  xs = [x/4.0 for x in mins]
  for alg in algs:
    sizes = [d[alg][avg][min][max][1].max/(avg * 1024.0) for min in mins]
    plt.plot(xs, sizes, label="alg=%s" % alg)
  ax = plt.gca()
  ax.set_ylim(bottom=0)
  ax.set_xlim(left=0, right=xs[-1])
  #plt.xscale('log')
  #plt.yscale('log')
  saveplt(FileName('sizemax','t', 'x', avg, max), 'size max vs min limit for avg=%sK,max=%sx' % (avg,max),
            'min limit fraction of avg', 'size max multiple of avg', xticks=xs)

def SizeDevVsMinLimit(data, avg, max):
  """Plot how size dev varys with min limit."""
  d = data
  xs = [x/4.0 for x in mins]
  for alg in algs:
    sizes = [d[alg][avg][min][max][1].dev/(avg * 1024.0) for min in mins]
    plt.plot(xs, sizes, label="alg=%s" % alg)
  ax = plt.gca()
  ax.set_ylim(bottom=0)
  ax.set_xlim(left=0, right=xs[-1])
  #plt.xscale('log')
  #plt.yscale('log')
  saveplt(FileName('sizedev', 't', 'x', avg, max) , 'size dev vs min limit for avg=%sK,max=%sx' % (avg,max),
            'min limit fraction of avg', 'size dev multiple of avg', xticks=xs)

def SizeAvgVsAvgTarget(data, alg, max):
  """Plot how size avg varys with avg target."""
  d = data
  xs = avgs
  for min in mins:
    sizes = [d[alg][avg][min][max][1].avg/(avg * 1024.0) for avg in avgs]
    plt.plot(xs, sizes, label="min=%sx" % (min/4.0))
  ax = plt.gca()
  #ax.set_ylim(bottom=0)
  #ax.set_xlim(left=0, right=xs[-1])
  plt.xscale('log')
  #plt.yscale('log')
  saveplt(FileName('sizeavg', alg, 't', 'x', max), 'size avg vs avg target for alg=%s,max=%sx' % (alg,max),
            'avg target KB', 'size avg multiple of avg', xticks=xs)

def PerfVsMinLimit(data, alg, max):
  """Plot how deduplication peformance varys with min limit."""
  d = data[alg]
  xs = [x/4.0 for x in mins]
  for avg in avgs:
    perfs = [d[avg][min][max][0]*100.0 for min in mins]
    plt.plot(xs, perfs, label="avg=%sK" % avg)
  ax = plt.gca()
  ax.set_ylim(bottom=0)
  ax.set_xlim(left=0, right=xs[-1])
  #plt.xscale('log')
  #plt.yscale('log')
  saveplt(FileName('perf', alg, 'x', 't', max), 'Deduplication vs min limit for alg=%s,max=%sx' % (alg,max),
            'min limit fraction of avg', 'found %', xticks=xs)

def PerfVsMaxLimit(data, alg, min):
  """Plot how deduplication peformance varys with max limit."""
  d = data[alg]
  xs = [float(x) for x in maxs]
  for avg in avgs:
    perfs = [d[avg][min][max][0]*100.0 for max in maxs]
    plt.plot(xs, perfs, label="avg=%sK" % avg)
  ax = plt.gca()
  ax.set_ylim(bottom=0)
  ax.set_xlim(left=0, right=xs[-1])
  #plt.xscale('log')
  #plt.yscale('log')
  saveplt(FileName('perf', alg, min, 't', 'x'), 'Deduplication vs max limit for alg=%s,min=%sx' % (alg,min),
            'max limit multiple of avg', 'found %', xticks=xs)

def PerfVsAvgSize(data, min, max):
  """Plot how deduplication peformance varys with avg target."""
  d = data
  xs = avgs
  for alg in algs:
    perfs = [d[alg][avg][min][max][0]*100.0 for avg in avgs]
    xs = [d[alg][avg][min][max][1].avg/1024.0 for avg in avgs]
    plt.plot(xs, perfs, label="alg=%s" % alg)
  ax = plt.gca()
  ax.set_ylim(bottom=0)
  ax.set_xlim(left=0, right=xs[-1])
  #plt.xscale('log')
  #plt.yscale('log')
  saveplt(FileName('perf', 't', min, 'x', max),
          'Deduplication vs avg target for min=%sx,max=%sx' % (min/4.0, max),
          'avg size KB', 'found %', xticks=avgs)

# This builds a results dict structured as;
# results[alg][avg][min][max] = (perf, blkstats, dupstats)
#algs = sorted(chunkers)
algs = ['weibull0', 'weibull1', 'weibull2', 'fastcdc']
data = {}
for alg in algs:
  data[alg] = GetFileData(alg)

alg0 = data[algs[0]]
avgs = sorted(alg0)
avg0 = alg0[avgs[0]]
mins = sorted(avg0)
min0 = avg0[mins[0]]
maxs = sorted(min0)

for alg in algs:
  PerfVsMinLimit(data, alg, maxs[-1])
  PerfVsMaxLimit(data, alg, mins[0])
  SizeAvgVsAvgTarget(data, alg, maxs[0])
  SizeAvgVsAvgTarget(data, alg, maxs[-1])
SizeMaxVsMinLimit(data, avgs[0], maxs[-1])
SizeDevVsMinLimit(data, avgs[0], maxs[-1])
PerfVsAvgSize(data, mins[0], maxs[-1])
PerfVsAvgSize(data, 0, 2)
PerfVsAvgSize(data, 0, 4)
PerfVsAvgSize(data, 1, 2)
PerfVsAvgSize(data, 1, 4)
PerfVsAvgSize(data, 1, 8)
PerfVsAvgSize(data, 2, 2)
PerfVsAvgSize(data, 2, 4)
PerfVsAvgSize(data, 2, 8)
PerfVsAvgSize(data, 3, 2)
