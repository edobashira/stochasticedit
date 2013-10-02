#! /usr/bin/env python

"""

Implementation sketch of Ristad & Yianilos' estimation procedure for
memoryless stochastic transducers.  See:

  Eric Sven Ristad and Peter N. Yianilos, `Learning String Edit Distance',
  IEEE Transactions on Pattern Analysis and Machine Intelligence 20(5):
  522--522, May 1998.

"""

__all__ = ('MemorylessStochasticTransducer',)


import math
import random
#import fst

try:
    import Numeric
    def _calloc_array(x, y):
        return Numeric.zeros((x, y), Numeric.Float)
except:
    def _calloc_array(x, y):
        return [[0.0 for j in xrange(y)] for i in xrange(x)]

_LOG_2 = math.log(2.0)


class MemorylessStochasticTransducer:

    """Class representing a stochastic two-tape transducer without
    any memory.  Defines a relation between strings over alphabet A
    and strings over alphabet B.  Stochastic edit distance between
    strings x and y is defined as the (negative log) probability
    assigned by the transducer to the pair (x,y).  Uses Ristad &
    Yianilos' method for estimating the (|A|+1)*(|B|+1) free
    parameters."""
    
    def __init__(self, A, B, epsilon, delta=None, stop=None, seed=1):
        """Creates a new transducer with alphabets A and B and empty
        string symbol epsilon.  Arc probabilities and the stopping
        probability may are optional parameters; if values are
        missing, random numbers are substituted and normalized to
        obtain a proper probability distribution."""
        # initialize the alphabets and the epsilon symbol
        self.A = A
        self.B = B
        self.epsilon = epsilon
        self.E = ([(a,b) for a in A for b in B]
                  + [(a,epsilon) for a in A]
                  + [(epsilon,b) for b in B])
        #print self.E
        # make sure we have a proper probability distribution
        self.R = random.Random(seed)
        self.delta = {}
        for e in self.E:
            if delta is None or not delta.has_key(e):
                self.delta[e] = self.R.random()
            else:
                self.delta[e] = float(delta[e])
        if stop is None:
            self.stop = 1#self.R.random()
        else:
            self.stop = float(stop)
        self.normalize()
        return

    def normalize(self):
        """Normalizes the parameters so that they sum to 1."""
        N = self.stop
        for e in self.E:
            N += self.delta[e]
        assert N > 0
        print "norms",self.stop,N
        self.stop /= N
        for e in self.E:
            self.delta[e] /= N
        return

    def stochastic_edit_distance(self, x, y):
        """Returns the negative log_2 probability of (x,y)."""
        alpha = self.forward_evaluate(x, y)
        return - math.log(alpha[len(x)][len(y)]) / _LOG_2

    def log_likelihood(self, C):
        """Returns the negative log likelihood of a collection of pairs."""
        ll = 0.0
        for x,y in C:
            ll += self.stochastic_edit_distance(x,y)
        return ll

    def forward_evaluate(self, x, y):
        """Returns a matrix of prefix probabilities."""
        T = len(x)
        V = len(y)
        alpha = _calloc_array(T+1, V+1)
        alpha[0][0] = 1.0
        for t in range(0, T+1):
            for v in range(0, V+1):
                # NB: The published paper says (v > 1) etc.
                # That would be wrong.
                if v > 0:
                    alpha[t][v] += (self.delta[(self.epsilon,y[v-1])]
                                    * alpha[t][v-1])
                if t > 0:
                    alpha[t][v] += (self.delta[(x[t-1],self.epsilon)]
                                    * alpha[t-1][v])
                if v > 0 and t > 0:
                    alpha[t][v] += (self.delta[(x[t-1],y[v-1])]
                                    * alpha[t-1][v-1])
        alpha[T][V] *= self.stop
        return alpha

    def backward_evaluate(self, x, y):
        """Returns a matrix of suffix probabilities."""
        T = len(x)
        V = len(y)
        beta = _calloc_array(T+1, V+1)
        beta[T][V] = self.stop
        for t in range(T, -1, -1):
            for v in range(V, -1, -1):
                if v < V:
                    beta[t][v] += (self.delta[(self.epsilon,y[v])]
                                   * beta[t][v+1])
                if t < T:
                    beta[t][v] += (self.delta[(x[t],self.epsilon)]
                                   * beta[t+1][v])
                if v < V and t < T:
                    beta[t][v] += (self.delta[(x[t],y[v])]
                                   * beta[t+1][v+1])
        return beta

    def expectation_maximization(self, C):
        """Performs a single EM iteration."""
        gamma_d = {}
        for z in self.E:
            gamma_d[z] = 0.0
        gamma = {'stop': 0.0, 'delta': gamma_d}
        for x,y in C:
            self.expectation_step(x, y, gamma)
        print "gammas",gamma
        self.maximization_step(gamma)
        print gamma
        return

    def expectation_step(self, x, y, gamma, scale=1.0):
        """Computes the expected counts of the edit operations generating
        (x,y) and adds them to gamma."""
        T = len(x)
        V = len(y)
        print "lengths",T,V
        alpha = self.forward_evaluate(x, y)
        if alpha[T][V] == 0.0:
            return
        beta = self.backward_evaluate(x, y)
        gamma['stop'] += scale
        gamma_d = gamma['delta']
        for t in range(0, T+1):
            for v in range(0, V+1):
                norm = scale * beta[t][v] / alpha[T][V]
                if t > 0:
                    z = (x[t-1], self.epsilon)
                    gamma_d[z] += alpha[t-1][v] * self.delta[z] * norm
                if v > 0:
                    z = (self.epsilon, y[v-1])
                    gamma_d[z] += alpha[t][v-1] * self.delta[z] * norm
                if t > 0 and v > 0:
                    z = (x[t-1], y[v-1])
                    gamma_d[z] += alpha[t-1][v-1] * self.delta[z] * norm
        return

    def maximization_step(self, gamma):
        """Sets the parameters to their expected counts and renormalizes."""
        # here's a neat way to update state in Python:
        self.__dict__.update(gamma)
        self.normalize()
        #print self.delta, self.stop
        return

def linear_chain(text, syms = None):
    """linear_chain(text, syms=None) -> linear chain acceptor for the given input text"""
    chain = fst.LogVectorFst();
    chain.start = chain.add_state()
    for i, c in enumerate(text):
        print i, i + 1, c
        chain.add_state()
        chain.add_arc(i, i + 1,  st[c], st[c], 0)
        #chain.add_arc(i, i+1, c)
    chain[i+1].final = True
    return chain

if __name__ == '__main__':
    phi = MemorylessStochasticTransducer('ht', 'ht', '<eps>')
    print phi.delta
    print phi.stop
    alpha = phi.forward_evaluate('hth', 'ttt')
    """
    st = fst.SymbolTable("syms")
    st["<eps>"] = 0
    st["h"] = 1
    st["t"] = 2
    f = fst.LogVectorFst()
    f.start = f.add_state()
    
    for a in phi.delta:
      i = st[a[0]]
      o = st[a[1]]
      w = -math.log(phi.delta[a])
      print i,o,w
      f.add_arc(0, 0, i, o, w)
    f[0].final = -math.log(phi.stop)

    f.write("E.ofst")
    X = linear_chain("hth", st)
    Y = linear_chain("ttt", st)
    XEY = X >> f >> Y
    XEY.write("XEY.ofst")"""
    beta = phi.backward_evaluate('hth', 'ttt')
    #assert math.fabs(alpha[3][2] - beta[0][0]) <= 1e-18
    print "alphas"
    print alpha
    print
    print "betas"
    print beta
    print
    print phi.delta, phi.stop
    print
    #exit(-1)
    C = [('hth','ttt')]
    for n in xrange(10):
        phi.expectation_maximization(C)
        print phi.log_likelihood(C)

    print phi.delta, phi.stop

## eof
