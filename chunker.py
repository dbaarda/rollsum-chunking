#!/usr/bin/pypy3 -O
"""
RollsumChunking modelling.

This will run tests for the specified chunker with different avg/min/max
length settings, and dump the summary data into a file in a directory.

Usage: %(cmd)s <chunker|weibull0|weibull1|weibull2|fastcdc|fastweibull2> [dir]
"""
from __future__ import print_function
import os
import pickle
import random
import sys
from math import e, gamma, log
from stats1 import Sample

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


class RandIter(object):
  """ A fast LCG random uint32 iterator.

  This also supports specifying a cycle period where it will repeat the
  values, and a fast skip(n) method for skipping over n values.
  """

  # Pythons random stuff is too slow, so we use a simple good-enough LCG
  # generator with modulus m = 2^32 for 32 bits. These values come from
  # Numerical Recipes.
  m = 2**32  # LCG modulus value.
  a = 1664525  # LCG multiplier value.
  c = 1013904223  # LCG increment value.
  b = m - 1  # fast bitmask version of m.

  def __init__(self, seed, cycle=2**32):
    self.seed=seed
    self.cycle=cycle
    self.value = seed
    self.count = 0  # count of total values produced.
    self.dat_n = 0  # count of values left this cycle.

  def __iter__(self):
    return self

  def __next__(self):
    if not self.dat_n:
      self.value = self.seed
      self.dat_n = self.cycle
    self.value = (self.a * self.value + self.c) & self.b
    self.count += 1
    self.dat_n -= 1
    return self.value

  # For python2 compatibility.
  next = __next__

  def skip(self, n):
    """ Skip over the next n random values. """
    self.count += n
    self.dat_n -= n
    # if past the cycle length, skip to the start of the last cycle.
    if self.dat_n < 0:
      self.value = self.seed
      self.dat_n %= self.cycle
      n = self.cycle - self.dat_n
    # https://www.nayuki.io/page/fast-skipping-in-a-linear-congruential-generator
    m, a, c = self.m, self.a, self.c
    a1 = self.a - 1
    ma = a1 * m
    self.value = (pow(a, n, m)*self.value + (pow(a, n, ma) - 1) // a1 * c) & self.b
    return self.value


class Data(object):
  """ Data source with rollsums and block hashes.

  It simulates a stream of data that starts with olen bytes of initial random
  data that is then repeated with modifications. The modifications are cycles
  of copied, inserted, and deleted data. The copy, insert, and delete have
  exponentially distributed random lengths with averages of clen, ilen, and
  dlen respectively.

  It simulates returning a 32bit rolling hash for each input byte with
  getroll(). A simulated strong hash of the previous block can be fetched with
  gethash(), which also starts a new block.

  It keeps counts of the number of bytes and duplicate bytes.
  """

  def __init__(self, olen, clen, ilen, dlen, seed=1):
    self.olen = olen
    self.clen = clen
    self.ilen = ilen
    self.dlen = dlen
    self.seed = seed
    # exponential distribution lambda parameters for clen/ilen/dlen.
    self.clambd = 1.0/clen
    self.ilambd = 1.0/ilen
    self.dlambd = 1.0/dlen
    self.reset()

  def reset(self):
    self.tot_c = 0         # total bytes scanned.
    self.dup_c = 0         # duplicate bytes scanned.
    self.blkh = 0          # the accumulated whole block hash.
    # Initialize the random generators for the original and inserted data.
    self.dat = RandIter(self.seed, self.olen)
    self.ins = RandIter(self.seed + 6)
    self.mod = random.Random(self.seed)
    self.cpystats = Sample()
    self.insstats = Sample()
    self.delstats = Sample()
    self.initcycle()

  def initcycle(self):
    self.cpy_n = int(self.mod.expovariate(self.clambd))
    self.ins_n = self.ilen and int(self.mod.expovariate(self.ilambd))
    self.del_n = self.dlen and int(self.mod.expovariate(self.dlambd))
    self.cpystats.add(self.cpy_n)
    self.insstats.add(self.ins_n)
    self.delstats.add(self.del_n)

  def getroll(self):
    if self.tot_c < self.olen:
      # Output initial data.
      h = self.dat.next()
    elif self.cpy_n:
      # Output copied data.
      h = self.dat.next()
      self.cpy_n -= 1
      self.dup_c += 1
    elif self.ins_n:
      # Output inserted data.
      h = self.ins.next()
      self.ins_n -= 1
    else:
      # do delete, setup next cycle, and recurse.
      self.dat.skip(self.del_n)
      self.initcycle()
      return self.getroll()
    # increment tot_c and update blkh.
    self.tot_c += 1
    self.blkh = hash((self.blkh, h))
    return h

  def gethash(self):
    """ Get a strong hash of the past l bytes and reset for a new block. """
    blkh, self.blkh = self.blkh, 0
    return blkh

  def __repr__(self):
    return "Data(olen=%s, clen=%s, ilen=%s, dlen=%s, seed=%s)" % (
        self.olen, self.clen, self.ilen, self.dlen, self.seed)

  def __str__(self):
    return "%r: tot=%d dup=%d(%4.1f%%)\n  cpy: %s\n  ins: %s\n  del: %s" % (
        self, self.tot_c, self.dup_c, 100.0 * self.dup_c / self.tot_c,
        self.cpystats, self.insstats, self.delstats)


class Chunker(object):
  """ A standard exponential chunker

  This is the standard simple chunker that gives an exponential distribution
  of block sizes between min and max. The only difference is it uses 'h<p`
  instead of 'h&mask==r' for the hash judgement, which supports arbitrary
  target block sizes, not just power-of-2 sizes. For tgt_len as the mean, the
  distribution's curves where x is measured from min_len and L is the
  normal exponential distribution lambda parameter are;

    f(x) = L
    CDF(x) = 1 - e^-(L*x)
    PDF(x) = L*e^-(L*x)
    mean = C + 1/L*(1-e^-(L*T))

  Where;

    L = 1/tgt_len
    C = min_len
    T = max_len - min_len

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
    if tgt_len <= 0:
      return min_len
    z = (max_len - min_len)/tgt_len
    return min_len + tgt_len * (1.0 - e**-z)

  @classmethod
  def get_tgt_len(cls, avg_len, min_len, max_len):
    """Get the tgt_len given an avg_len."""
    return solve(lambda x: cls.get_avg_len(x, min_len, max_len) - avg_len, 0.0)

  def reset(self):
    self.blocks = {}
    self.blkstats = Sample()
    self.dupstats = Sample()
    self.prob = 2**32 // self.tgt_len
    self.initblock()

  def initblock(self):
    self.blk_len = 0

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
    return "%s(tgt_len=%s, min_len=%s, max_len=%s)" % (
        self.__class__.__name__, self.tgt_len, self.min_len, self.max_len)

  def __str__(self):
    return "%r: avg_len=%s\n  blks: %s\n  dups: %s" % (
        self, self.avg_len, self.blkstats, self.dupstats)


class Weibull0Chunker(Chunker):
  """ Weibull0Chunker class.

  This uses a chunking probability criteria where the hash is treated as a
  fixed point number in the range 0.0 -> 1.0 and compared to a slowly
  increasing probability. The position x past min_len is a chunk boundary if h
  < f(x) where the f(x) "hazard function" is a function of x^P. This gives a
  Weibull block length distribution. For tgt_len as the mean and P curve
  power, the distribution's curves where x is measured from min_len and k and
  L are the normal Weibull parameters are;

    f(x) = M*x^P
    CDF(x) = 1 - e^-(M/k*x^k)
    PDF(x) = M*x^(k-1) * e^-(M/k*x^k)
    mean = C + L*gammalower((k+1)/k,(T/L)^k) + T*e^-((T/L)^k)

  Where;

    k = P + 1
    L = tgt_len/gamma(1+1/k)
    M = k/L^k = b*k
    C = min_len
    T = max_len - min_len

  The tgt_len for this chunker represents the distribution mean, not including
  the effects of min_len and max_len. This class uses P=0 (k=1) which makes it
  the same as a classic chunker, but subclasses can overide P for different
  variants.
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
    super(Weibull0Chunker, self).reset()
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

class Weibull1Chunker(Weibull0Chunker):
  P = 1

class Weibull2Chunker(Weibull0Chunker):
  P = 2


class WeibullT0Chunker(Weibull0Chunker):
  """WeibullT0 Chunker class.

  This is the similar to the Weibull Chunker except that min_len doesn't just
  shift the distribution to the right, instead it zero's the hazard function.
  This changes the distribution so it's not actually a Weibull distribution
  any more, unless min_len=0. The distribution's curves, where x is measured
  from min_len and k and L are the normal Weibull parameters, are;

    f(x) = M*(x+C)^P
    CDF(x) = 1 - e^-(M/k*((x+C)^k - C^k))
    PDF(x) = M*(x+C)^(k-1) * e^-(M/k*((x+C)^k - C^k))
    mean = L*e^((C/L)^k) * (gammalower((k+1)/k, ((T+C)/L)^k) -
        gammalower((k+1)/k, (C/L)^k)) + (C+T)*e^-(((T+C)/L)^k - (C/L)^k)

  Where;

    k = P + 1
    L = tgt_len/gamma(1+1/k)
    M = k/L^k
    C = min_len
    T = max_len - min_len

  The tgt_len for this chunker represents the weibull distribution mean, not
  including the effects of min_len and max_len. This class uses P=0 (k=1)
  which makes it the same as a classic chunker, but subclasses can overide P
  for different variants.
  """

  @classmethod
  def get_avg_len(cls, tgt_len, min_len, max_len):
    if tgt_len <= 0:
      return min_len
    k = cls.P + 1
    s = 1.0 + 1.0/k
    L = tgt_len / gamma(s)
    zc = (min_len/L)**k
    zt = (max_len/L)**k
    return L * e**zc * (gammalower(s, zt) - gammalower(s, zc)) + max_len*e**(zc-zt)

  def incblock(self):
    self.blk_len += 1
    x = self.blk_len
    if x >= self.min_len:
      self.prob = int(self.M * x**self.P)

class WeibullT1Chunker(WeibullT0Chunker):
  P = 1

class WeibullT2Chunker(WeibullT0Chunker):
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
  """ FastWeibull2Chunker class.

  This is the same as Weibull2Chunker except it aproximates it using integers only
  to increment prob every dx iterations. This turns out slower in Python, but
  it would almost certainly be a faster in C.
  """

  def reset(self):
    # Set incr needed by initblock() called by super().
    self.incr = 0
    super(FastWeibull2Chunker, self).reset()
    # set default update interval and scaling values. These scaling values ensure
    # that incr is large enough to be accurate enough (greater than 2^7).
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


def addtest(table, data, dsize, bsize, cls, bavg, bmin, bmax):
  try:
    table[bavg][bmin][bmax]
  except KeyError:
    bavg_len = bsize * bavg
    bmin_len = int(bavg_len * bmin)
    bmax_len = int(bavg_len * bmax)
    data.reset()
    chunker = cls.from_avg(bavg_len, bmin_len, bmax_len)
    result = runtest(chunker, data, 2*dsize)
    tableadd(table, result, bavg, bmin, bmax)


def alltests(cls, tsize, bsize):
  """Get results for different avg,min,max chunker args."""
  results = {}
  # Data size is tsize times the average 8*bsize blocks.
  dsize = tsize*bsize*8
  data = Data(dsize, bsize*16, bsize*8, bsize*4)
  for bavg in (1,2,4,8,16,32,64):
    for bmin in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7):
      addtest(results, data, dsize, bsize, cls, bavg, bmin, 8.0)
    for bmax in (1.5, 2.0, 4.0, 8.0):
      addtest(results, data, dsize, bsize, cls, bavg, 0.0, bmax)
    addtest(results, data, dsize, bsize, cls, bavg, 0.5, 2.0)
    addtest(results, data, dsize, bsize, cls, bavg, 0.2, 4.0)
    addtest(results, data, dsize, bsize, cls, bavg, 0.5, 4.0)
  return (tsize, bsize, results)


chunkers = dict(
    chunker=Chunker,
    weibull0=Weibull0Chunker,
    weibull1=Weibull1Chunker,
    weibull2=Weibull2Chunker,
    weibullt0=WeibullT0Chunker,
    weibullt1=WeibullT1Chunker,
    weibullt2=WeibullT2Chunker,
    fastcdc=FastCDCChunker,
    fastweibull2=FastWeibull2Chunker)

# Get the standard and (suspected) optimal chunker min fraction.
minchunkerstd = 0.2
minchunkeropt = log(2) / (1 + log(2))
# Get the fastcdc standard avg and min fraction (min=8K, tgt=4K+min, max=64k).
avgfastcdcstd = FastCDCChunker.get_avg_len(12*1024, 8*1024, 64*1024)
minfastcdcstd = 8.0*1024 / avgfastcdcstd
# Get the fastcdc optimized avg and min fraction (min=4K, tgt=4K+min, max=64k).
avgfastcdcopt = FastCDCChunker.get_avg_len(8*1024, 4*1024, 64*1024)
minfastcdcopt = 4.0*1024 / avgfastcdcstd


def usage(code, error=None, *args):
  if error:
    print(error % args)
  print(__doc__ % dict(cmd=os.path.basename(sys.argv[0])))
  sys.exit(code)


if __name__ == '__main__':
  cmd = sys.argv[1] if len(sys.argv) > 1 else None
  dir = sys.argv[2] if len(sys.argv) > 2 else '.'
  if cmd in ("-?", "-h", "--help", None):
    usage(0)
  if cmd not in chunkers:
    usage(1, "Error: invalid chunker argument %r.", cmd)
  cls = chunkers[cmd]
  results = alltests(cls, tsize=10000, bsize=1024)
  pickle.dump(results, open('%s/%s.dat' % (dir,cmd), 'wb'))
