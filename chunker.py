#!/usr/bin/pypy -O
"""
RollsumChunking modelling.

This will run tests for the specified chunker with different avg/min/max
length settings.

Usage: %(cmd)s [chunker|weibull0|weibull1|weibull2|fastcdc|fastweibull2]
"""
from __future__ import print_function
from math import *
from stats1 import Sample
import os
import random
import sys

def solve(f, x0=-1.0e9, x1=1.0e9, e=1.0e-9):
  """ Solve f(x)=0 for x where x0<=x<=x1 within +-e. """
  y0, y1 = f(x0), f(x1)
  # y0 and y1 must have different sign.
  assert y0*y1 <= 0
  while (x1 - x0) > e:
    xm = (x0 + x1) / 2.0
    ym = f(xm)
    if y0*ym > 0:
      x0,y0 = xm,ym
    else:
      x1,y1 = xm,ym
  return x0

def gammalower(s, z):
  # Make sure s is a float.
  s=float(s)
  # For large z it converges to gamma(s).
  if z >= 32:
    return gamma(s)
  tot = term = z**s * e**-z / s
  # For the range of z and s values we care about, this is enough iterations.
  for i in range(1,int(2*z)+12):
    term *= z / (s+i)
    tot += term
  return tot


class Data(object):
  """ Data source with rollsums and block hashes.

  It simulates a stream of data that starts with bnum blocks of initial random
  data that is then repeated every bnum blocks with modifications. The
  modifications consist of a modification every mnum blocks that starts 1/3 of
  a block into the first block, and replaces the next 1/7 of a block with 1/5
  of a block of new random data.

  It simulates returning a 32bit rolling hash for each input byte with
  getroll(). A simulated strong hash of the previous block can be fetched with
  gethash(), which also starts a new block.

  It keeps counts of the number of bytes and duplicate bytes.
  """

  def __init__(self, bsize=1024, bnum=512, mnum=5, seed=1):
    self.bsize = bsize     # block size.
    self.bnum = bnum       # number of blocks before repeating.
    self.mnum = mnum       # number of repeated blocks per change.
    self.seed = seed
    self.dat_p = bsize * bnum # period over which data repeats.
    self.mod_p = bsize * mnum # period over which changes happen.
    self.mod_o = bsize // 3    # offset at which changes happen.
    self.del_c = bsize // 7    # bytes deleted each change.
    self.ins_c = bsize // 5    # bytes inserted each change.
    self.mod_e = self.mod_o + self.ins_c
    self.reset()

  def reset(self):
    self.tot_c = 0         # total bytes scanned.
    self.dup_c = 0         # duplicate bytes scanned.
    self.blkh = 0          # the accumulated whole block hash.
    # Initialize the random generators for the original and inserted data.
    self.dat = random.Random(self.seed)
    self.ins = random.Random(self.seed + 6)

  def getroll(self):
    """ Get the next rolling hash. """
    c = self.tot_c
    # Start duplicating data every dat_p bytes.
    if (c % self.dat_p) == 0:
      self.dat.seed(self.seed)
    # After the first dat_p bytes, start making changes.
    if c < self.dat_p:
      h = self.dat.randrange(2**32)
    else:
      # Set i to the offset past the periodic modify point.
      i = c % self.mod_p
      # At offset mod_o modify stuff.
      if i == self.mod_o:
        # delete del_c bytes by sucking them out of dat.
        for d in range(self.del_c):
          self.dat.randrange(2**32)
        #print "%12d: start replace, del=%d" % (self.tot_c, self.del_c)
      #elif i == self.mod_e:
      #  print "%12d: stop replace, ins=%d" % (self.tot_c, self.ins_c)
      # Between mod_o and mod_e insert new data, otherwise use duplicate data.
      if self.mod_o <= i < self.mod_e:
        h = self.ins.randrange(2**32)
      else:
        self.dup_c += 1
        h = self.dat.randrange(2**32)
    self.tot_c += 1
    # update blkh.
    self.blkh = hash((self.blkh, h))
    return h

  def gethash(self):
    """ Get a strong hash of the past l bytes and reset for a new block. """
    blkh, self.blkh = self.blkh, 0
    return blkh

  def __repr__(self):
    return "Data(bsize=%s, bnum=%s, mnum=%s, seed=%s)" % (self.bsize, self.bnum, self.mnum, self.seed)

  def __str__(self):
    return "%r: tot=%d dup=%d(%4.1f%%)" % (
        self, self.tot_c, self.dup_c, 100.0 * self.dup_c / self.tot_c)


class Chunker(object):
  """ A standard exponential chunker

  This is the standard simple chunker that gives an exponential distribution
  of block sizes between min and max. The only difference is it uses 'h<p`
  instead of 'h&mask==r' for the hash judgement, which supports arbitrary
  target block sizes, not just power-of-2 sizes.

  The tgt_len for this chunker represents the exponential distribution mean
  size, not including the affects of min_len and max_len.
  """

  MIN_LEN, MAX_LEN = 0, 2**32

  def __init__(self, tgt_len, min_len=MIN_LEN, max_len=MAX_LEN):
    assert min_len < max_len
    self.tgt_len = tgt_len
    self.min_len = min_len
    self.max_len = max_len
    self.avg_len = self.get_avg_len(tgt_len, min_len, max_len)
    self.reset()

  @classmethod
  def from_avg(cls, avg_len, min_len=MIN_LEN, max_len=MAX_LEN):
    """Initialize using the avg_len."""
    tgt_len = int(cls.get_tgt_len(avg_len, min_len, max_len) + 0.5)
    return cls(tgt_len, min_len, max_len)

  @classmethod
  def get_avg_len(cls, tgt_len, min_len, max_len):
    """Get the avg_len given a tgt_len."""
    try:
      return min_len + (1.0 - e**(float(min_len - max_len)/tgt_len)) * tgt_len
    except ZeroDivisionError:
      return min_len

  @classmethod
  def get_tgt_len(cls, avg_len, min_len, max_len):
    """Get the tgt_len given an avg_len."""
    return solve(lambda x: cls.get_avg_len(x, min_len, max_len) - avg_len, 0.0)

  def reset(self):
    self.blocks = {}
    self.blkstats = Sample()
    self.dupstats = Sample()
    self.initblock()

  def initblock(self):
    self.blk_len = 0
    self.prob = 2**32 // self.tgt_len

  def incblock(self):
    self.blk_len += 1

  def isblock(self, r):
    """ Checks if rollsum r is a break point and increments the block. """
    self.incblock()
    return self.blk_len >= self.min_len and (r < self.prob or self.blk_len >= self.max_len)

  def addblock(self, h):
    """ Adds a block with hash h and initializes for the next block. """
    l = self.blk_len
    b = (h, l)
    n = self.blocks[b] = self.blocks.setdefault(b, 0) + 1
    self.blkstats.add(l)
    if n > 1:
      self.dupstats.add(l)
    self.initblock()

  def __repr__(self):
    return "%s(tgt_len=%s, min_len=%s, max_len=%s)" % (self.__class__.__name__, self.tgt_len, self.min_len, self.max_len)

  def __str__(self):
    return "%r: avg_len=%s\n  blks: %s\n  dups: %s" % (self, self.avg_len, self.blkstats, self.dupstats)


class WeibullChunker(Chunker):
  """ WeibullChunker class.

  This uses a chunking probability criteria where the hash is treated as a
  fixed point number in the range 0.0 -> 1.0 and compared to a slowly
  increasing probability. The position x is a chunk boundary if h < p where
  the p "hazard function" is;

  p = M * x^P

  This gives a Weibull block length distribution. For tgt_len as the mean and
  P curve power, the Weibull k and L (lambda) parameters, and the resulting M
  values are;

  k = P + 1
  L = tgt_len/gamma(1+1/k)
  M = k/L^k = b*k

  The tgt_len for this chunker represents the distribution mean, not
  including the effects of min_len and max_len.

  This class uses P=0 which makes it the same as a classic chunker, but
  subclasses can overide P for different variants.
  """
  P = 0

  @classmethod
  def get_avg_len(cls, tgt_len, min_len, max_len):
    # Getting the average length for the distribution chopped off at 't' is
    # calculated as follows;
    #
    # avg = integ(x*PDF(x), 0, t) + t*(1-CDF(t))
    # x*PDF(x) = k*(x/L)^k*e^(-(x/L)^k)
    # t*(1-CDF(t)) = t*e^(-(t/L)^k)
    # avg = integ(k*(x/L)^k*e^(-(x/L)^k),0,t) + t*e^(-(t/L)^k)
    #     = L*gammalower(1 + 1/k, (t/L)^k) + t*e^(-(t/L)^k)
    if tgt_len <= 0:
      return min_len
    k = cls.P + 1
    s = 1.0 + 1.0/k
    L = tgt_len / gamma(s)
    t = max_len - min_len
    z = (t/L)**k
    return min_len + L * gammalower(s, z) + t*e**(-z)

  def reset(self):
    super(WeibullChunker, self).reset()
    # Set the M probability multiplier.
    k = self.P + 1
    L = self.tgt_len / gamma(1.0 + 1.0/k)
    self.M = 2**32 * k / L**k

  def initblock(self):
    self.blk_len = 0
    self.prob = 0.0

  def incblock(self):
    self.blk_len += 1
    x = self.blk_len - self.min_len
    if x > 0:
      self.prob = int(self.M * x**self.P)

class Weibull1Chunker(WeibullChunker):
  P = 1

class Weibull2Chunker(WeibullChunker):
  P = 2


class FastCDCChunker(Chunker):
  """ FastCDCChunker class.

  This implements FastCDC's chunking algorithm modified to use uses a 'h<p'
  hash judgment to support arbitrary tgt_len values.

  The tgt_len for this chunker is the length where the probability steps up
  from 1/4x to 4x the normal exponential distribution probablity. Note that
  min_len doesn't offset this.
  """

  @classmethod
  def get_avg_len(cls, tgt_len, min_len, max_len):
    if tgt_len <= min_len:
      # It's the same as Chunker with tgt_len/4.
      return Chunker.get_avg_len(tgt_len // 4, min_len, max_len)
    if tgt_len >= max_len:
      # It's the same as Chunker with tgt_len*4.
      return Chunker.get_avg_len(tgt_len * 4, min_len, max_len)
    tgt_s = tgt_len * 4
    tgt_l = tgt_len // 4
    s_avg = min_len + tgt_s # avg length of first exponential distribution.
    l_avg = tgt_len + tgt_l # avg length of second exponential distribution.
    s_cut_avg = tgt_len + tgt_s # avg length of first dist cut after tgt_len.
    l_cut_avg = max_len + tgt_l # avg length of second dist cut after max_len.
    s_cut_fract = e**(float(min_len - tgt_len)/tgt_s) # fraction of first dist cut.
    l_cut_fract = e**(float(tgt_len - max_len)/tgt_l) # fraction of second dist cut.
    return s_avg + s_cut_fract*(l_avg - s_cut_avg + l_cut_fract*(max_len - l_cut_avg))

  def initblock(self):
    self.blk_len = 0
    # prob is 1/4 of the normal tgt_len probability.
    self.prob = 2**30 // self.tgt_len

  def incblock(self):
    self.blk_len += 1
    if self.blk_len == self.tgt_len:
      # prob is 4x the normal tgt_len probability.
      self.prob = 2**34 // self.tgt_len


class FastWeibull2Chunker(Weibull2Chunker):
  """ FastWeibullChunker class.

  This is the same as Weibull2Chunker except it aproximates it using integers only
  to increment prob every dx iterations. This turns out slower in Python, but
  it would almost certainly be a faster in C.
  """

  def reset(self):
    # Set incr needed by initblock() called by super().
    self.incr = 0
    super(FastWeibull2Chunker, self).reset()
    # set default update interval and scaling values. These scaling values ensure
    # that incr is large enough to be accurate (greater than 2^7).
    dx, incr = 1, 2*self.M
    while incr < 128:
      dx *= 2
      incr = 2 * self.M * dx**2
    self.incr, self.mask = int(incr + 0.5), dx - 1
    # Call initblock() to initialize again with correctly set incr.
    self.initblock()

  def initblock(self):
    self.blk_len = 0
    self.prob = 0
    self.step = self.incr // 2

  def incblock(self):
    self.blk_len += 1
    x = self.blk_len - self.min_len
    if (x > 0) and (x & self.mask == 0):
      self.prob += self.step
      self.step += self.incr


def runtest(chunker, data, data_len):
  blocks = {}
  # Stop after we've read enough data and finished a whole block.
  while data.tot_c < data_len or chunker.blk_len:
    if chunker.isblock(data.getroll()):
      chunker.addblock(data.gethash())
  print(data)
  print(chunker)
  assert data.tot_c == chunker.blkstats.sum
  tot_n, dup_n = chunker.blkstats.num, chunker.dupstats.num
  tot_c, dup_c = chunker.blkstats.sum, chunker.dupstats.sum
  perf = float(dup_c) / data.dup_c
  print("bytes: tot=%s dup=%s(%4.2f%%)" % ( tot_c, dup_c, 100.0 * dup_c / tot_c))
  print("blocks: tot=%s dup=%s(%4.2f%%)" % ( tot_n, dup_n, 100.0 * dup_n / tot_n))
  print("found: %4.2f%%" % (100.0 * perf))
  print()
  return perf, chunker.blkstats, chunker.dupstats


def tableadd(table, value, *args):
  # Adds an entry to a nested dict of dicts keyed by the *args.
  for k in args[0:-1]:
    table = table.setdefault(k, {})
  table[args[-1]] = value


def alltests(cls, tsize, bsize):
  """Get results for different avg,min,max chunker args."""
  results = {}
  data = Data(bnum=tsize, bsize=bsize, mnum=4)
  for bavg in (1,2,4,8,16,32,64):
    bavg_len = bavg * 1024
    for bmin in (0, 1, 2, 3):
      bmin_len = bavg_len * bmin // 4
      for bmax in (16, 8, 4, 2):
        bmax_len = bavg_len * bmax
        data.reset()
        chunker = cls.from_avg(bavg_len, bmin_len, bmax_len)
        result = runtest(chunker, data, 2*tsize*bsize)
        tableadd(results, result, bavg, bmin, bmax)
  return results


chunkers = dict(
    chunker=Chunker,
    weibull0=WeibullChunker,
    weibull1=Weibull1Chunker,
    weibull2=Weibull2Chunker,
    fastcdc=FastCDCChunker,
    fastweibull2=FastWeibull2Chunker)

def usage(code, error=None, *args):
  if error:
    print(error % args)
  print(__doc__ % dict(cmd=os.path.basename(sys.argv[0])))
  sys.exit(code)

if __name__ == '__main__':
  cmd = sys.argv[1] if len(sys.argv) > 1 else None
  if cmd in ("-?", "-h", "--help", None):
    usage(0)
  if cmd not in chunkers:
    usage(1, "Error: invalid chunker argument %r.", cmd)
  cls = chunkers[cmd]
  alltests(cls, tsize=1000, bsize=8*1000)
