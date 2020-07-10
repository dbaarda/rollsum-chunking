#!/usr/bin/pypy
"""
RollsumChunking modelling.

This uses a chunking probability criteria where the hash is treated as a fixed
point number in the range 0.0 -> 1.0 and compared to a probablity of;

p = 2*(i - min_size)^2 / target_size^3

The position i is a chunk boundary if h <= p.
"""
from stats1 import Sample
import random


class Data(object):
  """ Data source with rollsums and block hashes. """

  def __init__(self, bsize=1024, bnum=512, mnum=5, seed=1):
    self.bsize = bsize     # block size.
    self.bnum = bnum       # number of blocks per cycle.
    self.mnum = mnum       # number of blocks per change.
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
    # Initialize the random generators.
    self.dat = random.Random(self.seed)
    self.ins = random.Random(self.seed +6)

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
        # delete del_c bytes...
        for d in range(self.del_c):
          self.dat.randrange(2**32)
      # Between mod_o and mod_e insert new data, otherwise use duplicate data.
      if self.mod_o <= i < self.mod_e:
        h = self.ins.randrange(2**32)
      else:
        self.dup_c += 1
        h = self.dat.randrange(2**32)
    self.tot_c += 1
    # update blkh.
    self.blkh += h
    return h

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

  def __init__(self, avg_len, min_len=None, max_len=None):
    if min_len == None:
      min_len = avg_len / 5
    tgt_len = avg_len - min_len
    if max_len == None:
      max_len = tgt_len * 4 + min_len
    assert min_len < avg_len < max_len
    self.avg_len = avg_len
    self.min_len = min_len
    self.max_len = max_len
    self.tgt_len = tgt_len
    self.reset()
    
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
    return "%s(avg_len=%s, min_len=%s, max_len=%s)" % (self.__class__.__name__, self.avg_len, self.min_len, self.max_len)

  def __str__(self):
    
    return "%r: %s" % (self, str(self.stats)[10:])


class NormChunker(Chunker):

  def reset(self):
    # step and incr are an incremental way to calculate p = K * x^2. The initial
    # step for incrementing p is p calculated for x=1, and at each update we
    # increment step by 2x that much.
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
      self.prob += self.step
      self.step += self.incr


class FastNormChunker(Chunker):

  def reset(self):
    # set default update interval and scaling values. These scaling values ensure
    # that incr is large enough to be accurate (greater than 2^7).
    self.di = 1  # update interval.
    self.incr = (2**32 * 4) / (self.tgt_len**3)
    while incr < 128:
      self.di *= 2
      self.incr = (2**32 * 4 * self.di**2) / (self.tgt_len**3)
    super(FastNormChunker, self).reset()
    
  def initblock(self):
    self.blk_len = 0
    self.prob = 0
    self.step = self.incr / 2

  def incblock(self):
    """ Returns the block length or -1 if not a break point. """
    self.blk_len += 1
    x = self.blk_len - self.min_len
    if (x > 0) and (x % self.di == 0):
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

tsize=2*1024
bsize=8*1024

data = Data(bnum=tsize/2, bsize=bsize, mnum=2)
for bavg in (1,2,4,8,16,32,64):
  bavg *= 1024
  for bmin in (0, bavg / 4, bavg / 2, bavg * 3 / 4):
    btgt = bavg - bmin
    for bmax in (16, 8, 4, 2):
      bmax = bmax*btgt + bmin
      data.reset()
      chunker = Chunker(bavg, bmin, bmax)
      runtest(data,chunker, tsize*bsize)
