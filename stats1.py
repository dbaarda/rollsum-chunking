#!/usr/bin/python
# Data structures and methods for analysis of data.

class Sample(object):
  
  def __init__(self):
    self.num = 0
    self.sum = 0.0
    self.sum2 = 0.0
    self.minv = float('inf')
    self.maxv = -float('inf')
  
  def add(self, v):
    self.num += 1
    self.sum += v
    self.sum2 += v*v
    if v > self.maxv:
      self.maxv = v
    if v < self.minv:
      self.minv = v
      
  @property
  def avg(self):
    return self.sum/self.num
  
  @property
  def var(self):
    return (self.sum2-(self.sum*self.sum/self.num))/(self.num)
  
  @property
  def dev(self):
    return self.var ** 0.5
  
  def __repr__(self):
    return "Sample()"

  def __str__(self):
    if self.num:
      return "%r: num=%d avg=%f dev=%f min=%s max=%s" % (
          self, self.num, self.avg, self.dev, self.minv, self.maxv)
    else:
      return "%r: num=0"

#class Histogram:
#       def Add(i):
#       def Inc(i,v):
#       def Count():
#       def High():
#       def Low():
#       def Best():
#       def Worst():
#       def Median():

if __name__ == '__main__':
  s = Sample()
  s.add(1.0)
  s.add(2.0)
  print s
