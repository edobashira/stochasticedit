#!/usr/bin/env python
import sys, random

class LogSemiring:
  def plus(w1, w2):
    return w1;

  def times(w1, w2):
    return w1

  def divide(w1, w2):
    return w1

  def one():
    return 0

  def zero():
    return inf

class RealSemiring:
  def plus(w1, w2):
    return w1 + w2;

  def times(w1, w2):
    return w1 * w1

  def divide(w1, w2):
    return w1 / w2

  def one():
    return 1

  def zero():
    return 0

class TropicalSemiring:
  def plus(w1, w2):
    return min(w1, w2)

  def times(w1, w2):
    return w1 + w2

  def divide(w1, w2):
    return w1 - w2

  def one():
    return 0

  def zero():
    return inf



class SymbolTable:
  def __init__(self):
    self.str2int = {}
    self.int2str = {}

  def num_symbols(self):
    return len(self.str2int)

  def find(self, a):
    if a in self.str2int:
      return self.str2int[a]
    return -1

  def find_sym(self, a):
    if a in self.int2str:
      return self. int2str[a]
    return "None"

  def add(self, a):
    if a not in self.str2int:
      sym = len(self.str2int)
      self.str2int[a] = sym
      self.int2str[sym] = a
      return sym
      #print >> sys.stderr, "adding ", a, sym
    else:
      return self.find(a)

  def dump(self):
    for k,v in self.str2int.items():
      print k,v


def initprobs(st):
  R = random.Random()
  S = st.num_symbols()
  C = [[ 0.1 for x in xrange(S)] for x in xrange(S)]
  e = 0
  H = st.find('H')
  T = st.find('T')
  
  #C[H][H] = 0.8
  #C[H][T] = 0.3
  #C[T][H] = 0.9
  #C[T][T] = 0.1
  
  #C[H][e] = 1.1
  #C[T][e] = 0.4
  
  #C[e][H] = 0.6
  #C[e][T] = 0.5
  
  #C[e][e] = 1

  """
  C[T][e] = 0.08345724793134081
  C[e][H] = 0.12098161795585818 
  C[H][T] = 0.15734347805599638 
  C[H][H] = 0.02494748152394602 
  C[e][T] = 0.14644268912418149 
  C[T][H] = 0.14181044459594713 
  C[H][e] = 0.091987699281077
  C[T][T] = 0.04735880329622476"""
  D = 0.1 #R.random()
  return C,D


def alphas(X, Y, C, D):
  #print "X,Y",X,Y
  M = len(X)
  N = len(Y)
  X = [0] + X
  Y = [0] + Y
  alpha = [[0 for x in xrange(N + 1)] for x in xrange(M + 1)] 
  alpha[0][0] = 1
  for m in range(1, M + 1):
    alpha[m][0] = C[X[m]][0] * alpha[m - 1][0] 
  for n in range(1, N + 1):
    alpha[0][n] = C[0][Y[n]] * alpha[0][n - 1] 
  for  m in range(1, M + 1):
    for n in range(1, N + 1):
      s = alpha[m - 1][n - 1] * C[X[m]][Y[n]]  
      d = alpha[m - 1][n]     * C[X[m]][0]
      i = alpha[m][n - 1]     * C[0][Y[n]]
      alpha[m][n] = s + d + i
  alpha[m][n] *= D
  return alpha


#input, output, alpha, beta, prob dist, num symbols
def EM(X, Y, A, B, C, S, lamba = 1):
  M = len(X)
  N = len(Y)
  X = [0] + X
  Y = [0] + Y
  #sigma = [[0 for x in xrange(N + 1)] for x in xrange(M + 1)] 
  A_m_n = A[-1][-1]
  #print A_m_n
  
  # Table of posteriors would be bigger because it is 
  # a grid with diagonal
  sigma = [[0 for x in xrange(S)] for x in xrange(S)]
  total = lamba
  final = 1
  #for m in range(1, M + 1):
    #x = X[m] #x symbol
    #print sigma[m][0]
    #print m,0,A[m - 1][0],C[x][0],B[m][0]
  
  #for n in range(1, N + 1):
    #y = Y[n] #y symbol
  #  sigma[0][n] 
  
  for m in range(0, M + 1):
    for n in range(0, N + 1):
      x = X[m] #x symbol
      y = Y[n] #x symbol
      if m > 0:
        p = A[m - 1][n] * C[x][0] * B[m][n] / A_m_n
        sigma[x][0] += p
        total += p
      if n > 0:
        p = A[m][n - 1] * C[0][y] * B[m][n] / A_m_n
        sigma[0][y] += p
        total += p
      if m > 0 and n > 0:
        p = A[m - 1][n - 1] * C[x][y] * B[m][n] / A_m_n
        sigma[x][y] += p
        total += p
  print "total",total
  for m in range(S):
    for n in range(S):
      print sigma[m][n]
      sigma[m][n] /= total
  final /= total
  #print final
  return sigma,final
        

def betas(X, Y, C, D):
  M = len(X)
  N = len(Y)
  X = [0] + X
  Y = [0] + Y
  beta = [[0 for x in xrange(N + 1)] for x in xrange(M + 1)] 
  beta[M][N] = D
  for m in range(M - 1, -1, -1):
    beta[m][N] = C[X[m + 1]][0] * beta[m + 1][N] 
  for n in range(N - 1, -1, -1):
    beta[M][n] = C[0][Y[n + 1]] * beta[M][n + 1] 
  for  m in range(M - 1, -1, -1):
    for n in range(N - 1, -1, -1):
      s = beta[m + 1][n + 1] * C[X[m + 1]][Y[n + 1]]  
      d = beta[m + 1][n]    * C[X[m + 1]][0]
      i = beta[m][n + 1]    * C[0][Y[n + 1]]
      beta[m][n] = s + d + i
  return beta


def print_matrix(M):
  for a in M:
    print " ".join("%.9f" % x for x in a)

def alloc_matrix(ist, ost):
  ni = ist.num_symbols() - 1
  ni = ost.num_symbols() - 1
  C = [[1.0 for x in xrange(ni)] for x in xrange(no)]
  return C

def print_probs(C, D, st):
  for c in xrange(len(C)):
    for d in xrange(len(C[c])):
      print "%s \t\t %s \t\t %.5f" % (st.find_sym(c), st.find_sym(d),C[c][d])
  print "final %.5f" % D
  print C,D


def main():
  X = "H".split(' ')
  Y = "H".split(' ')
  st = SymbolTable()
  st.add("<eps>")
  #print X,Y
  
  for i in range(len(X)):
    X[i] = st.add(X[i])
  for i in range(len(Y)):
    Y[i] = st.add(Y[i])
  #st.dump()
  print st.num_symbols(),"symbols"
  C, D = initprobs(st)
  print_probs(C,D,st)
  alpha = alphas(X, Y, C, D)
  #print X,Y
  #print alpha
  #beta = betas(X, Y, C, D)
  #print 
  #print_matrix(alpha)
  #print
  #print_matrix(beta)
  #print
  for i in range(0, 1):
    alpha = alphas(X, Y, C, D)
    print_matrix(alpha)
    beta = betas(X, Y, C, D)
    C, D = EM(X, Y, alpha, beta, C, 2)
    C[0][0] = 1
    print_probs(C,D,st)
      # print alpha[-1][-1]
  #print beta[-1][-1]

if __name__ == "__main__":
  main()
