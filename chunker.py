#!/usr/bin/pypy -O
"""
RollsumChunking modelling.

"""
from math import *
from stats1 import Sample
import random

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

  It keeps counts of the number of bytes, duplicate bytes, and blocks fetched.
  """

  def __init__(self, bsize=1024, bnum=512, mnum=5, seed=1):
    self.bsize = bsize     # block size.
    self.bnum = bnum       # number of blocks before repeating.
    self.mnum = mnum       # number of repeated blocks per change.
    self.seed = seed
    self.dat_p = bsize * bnum # period over which data repeats.
    self.mod_p = bsize * mnum # period over which changes happen.
    self.mod_o = bsize / 3    # offset at which changes happen.
    self.del_c = bsize / 7    # bytes deleted each change.
    self.ins_c = bsize / 5    # bytes inserted each change.
    self.mod_e = self.mod_o + self.ins_c
    self.reset()

  def reset(self):
    self.tot_c = 0         # total bytes scanned.
    self.dup_c = 0         # duplicate bytes scanned.
    self.blk_n = 0         # number of block hashes returned.
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
    self.addhash(h)
    return h

  def addhash(self, h):
    self.blkh = hash((self.blkh, h))

  def gethash(self, l):
    """ Get a strong hash of the past l bytes and reset for a new block. """
    blkh, self.blkh = self.blkh, 0
    self.blk_n += 1
    return blkh

  def __repr__(self):
    return "Data(bsize=%s, bnum=%s, mnum=%s, seed=%s)" % (self.bsize, self.bnum, self.mnum, self.seed)

  def __str__(self):
    return "%r: tot=%d dup=%d(%4.1f%%), blk=%d" % (
        self, self.tot_c, self.dup_c, 100.0 * self.dup_c / self.tot_c, self.blk_n)


class Chunker(object):
  """ A standard exponential chunker

  This is the standard simple chunker that gives an exponential distribution
  of block sizes between min and max. The only difference is it uses 'h<=p`
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
    tgt_len = cls.get_tgt_len(avg_len, min_len, max_len)
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
    return solve(lambda x: cls.get_avg_len(x, min_len, max_len) - avg_len, min_len, max_len)

  def reset(self):
    self.stats = Sample()
    self.initblock()

  def initblock(self):
    self.blk_len = 0
    self.prob = 2**32 / self.tgt_len

  def incblock(self):
    self.blk_len += 1

  def getblock(self, h):
    """ Returns the block length or 0 if not a break point. """
    self.incblock()
    if self.blk_len >= self.min_len and (h < self.prob or self.blk_len == self.max_len):
      l = self.blk_len
      self.stats.add(l)
      self.initblock()
      return l
    return 0

  def __repr__(self):
    return "%s(tgt_len=%s, min_len=%s, max_len=%s)" % (self.__class__.__name__, self.tgt_len, self.min_len, self.max_len)

  def __str__(self):
    return "%r: avg_len=%s %s" % (self, self.avg_len, self.stats)


class NormChunker(Chunker):
  """ NormChunker class.

  This uses a chunking probability criteria where the hash is treated as a fixed
  point number in the range 0.0 -> 1.0 and compared to a probablity of;

  p = 2*(i - min_size)^2 / target_size^3

  The position i is a chunk boundary if h < p.

  The tgt_len for this chunker represents the distribution mode point, not
  including the effects of min_len and max_len.
  """

  @classmethod
  def get_avg_len(cls, tgt_len, min_len, max_len):
    assert min_len + 2*tgt_len <= max_len
    return min_len + tgt_len * (3.0/2.0)**(1.0/3.0) * gamma(4.0/3.0)

  @classmethod
  def get_tgt_len(cls, avg_len, min_len, max_len):
    assert 2*avg_len - min_len <= max_len
    return (avg_len - min_len) * ((2.0/3.0)**(1.0/3.0)) / gamma(4.0/3.0)

  def reset(self):
    # self.K = 2**32 * 2.0 / self.tgt_len**3
    # step and incr are an incremental way to calculate p = K * x^2. The initial
    # step for incrementing p is p calculated for x=1, and at each update we
    # increment step by 2x that much. We calculate the increment here since it
    # is fixed, and reset step to half that in initblock().
    self.incr = 2**32 * 4.0 / self.tgt_len**3
    super(NormChunker, self).reset()

  def initblock(self):
    self.blk_len = 0
    self.prob = 0.0
    self.step = self.incr / 2.0

  def incblock(self):
    """ Returns the block length or -1 if not a break point. """
    self.blk_len += 1
    if self.blk_len > self.min_len:
      # self.prob = self.K * (self.blk_len - self.min_len)**2
      self.prob += self.step
      self.step += self.incr


class FastNormChunker(NormChunker):
  """ FastNormChunker class.

  This is the same as NormChunker except it aproximates it using integers only
  to increment prob every dx iterations. This turns out slower in Python, but
  it would almost certainly be a faster in C.
  """

  def reset(self):
    # set default update interval and scaling values. These scaling values ensure
    # that incr is large enough to be accurate (greater than 2^7).
    dx, incr = 1, (2**32 * 4) / (self.tgt_len**3)
    while incr < 128:
      dx *= 2
      incr = (2**32 * 4 * dx**2) / (self.tgt_len**3)
    self.incr, self.mask = incr, dx - 1
    super(FastNormChunker, self).reset()

  def initblock(self):
    self.blk_len = 0
    self.prob = 0
    self.step = self.incr / 2

  def incblock(self):
    """ Returns the block length or 0 if not a break point. """
    self.blk_len += 1
    x = self.blk_len - self.min_len
    if (x > 0) and (x & self.mask == 0):
      self.prob += self.step
      self.step += self.incr


def runtest(data, chunker, data_len):
  blocks = {}
  l = 0
  # Stop after we've read enough data and finished a whole block.
  while data.tot_c < data_len or l == 0:
    h = data.getroll()
    l = chunker.getblock(h)
    if l:
      b = (data.gethash(l), l)
      blocks[b] = blocks.setdefault(b, 0) + 1
      #print "%12d: %08x %016x len=%d dup=%s" % (data.tot_c, h, b[0], b[1], blocks[b]-1)
  # get stats on duplicate blocs.
  tot_c = tot_n = 0
  dup_c = dup_n = 0
  for (h,l), n in blocks.items():
    tot_n += n
    tot_c += n * l
    dup_n += (n - 1)
    dup_c += (n - 1) * l
  print data
  print chunker
  print "bytes: tot=%s dup=%s(%4.2f%%)" % ( tot_c, dup_c, 100.0 * dup_c / tot_c)
  print "blocks: tot=%s dup=%s(%4.2f%%)" % ( tot_n, dup_n, 100.0 * dup_n / tot_n)
  print "found: %4.2f%%" % (100.0 * dup_c / data.dup_c)
  print
  return tot_n, tot_c, dup_n, dup_c


tsize=2*1000
bsize=8*1000
data = Data(bnum=tsize/2, bsize=bsize, mnum=4)
for bavg in (1,2,4,8,16,32,64):
  bavg *= 1024
  for bmin in (0, bavg / 4, bavg / 2, bavg * 3 / 4):
    for bmax in (16, 8, 4, 2):
      bmax = bmax*bavg
      data.reset()
      chunker = NormChunker.from_avg(bavg, bmin, bmax)
      runtest(data, chunker, tsize*bsize)

# Dimensions for graphs;
#
# chunker (chunker, normchunker)
# min size (0, 1/4, 1/2, 3/4)
# max_size (16, 8, 4, 2)
# avg_size (1,2,4, 8, 16, 32, 64)
#
#
