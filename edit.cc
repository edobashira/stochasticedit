/*
 * Implementation of 
 * 
 *
 *
 */

#include <iostream>
#include <fst/vector-fst.h>
#include <fst/compose.h>
#include <fst/shortest-distance.h>

//#include "thread_pool.h"
using namespace std;
using namespace fst;

typedef VectorFst<LogArc> LogVectorFst;

template<class Arc>
void InitEditFst(VectorFst<Arc>* fst, int numsyms) {
  typedef typename Arc::Weight W;
  typedef typename Arc::Label L;
  typedef typename Arc::StateId S;

  //[eps, H, T] = [0, 1, 2]
  float probs[3][3] = { {1,   0.6, 0.5},
                        {1.1, 0.8, 0.3},
                        {0.4, 0.9, 0.1} };
  const char *strings[] = {"<eps>","H","T"};
  fst->DeleteStates();
  int st = fst->AddState();
  fst->SetStart(st);
  fst->SetFinal(st, W::One());
  for (int i = 1; i < 3; ++i) {
    fst->AddArc(st, Arc(i, 0, log(probs[i][0]), st));
    fst->AddArc(st, Arc(0, i, log(probs[0][i]), st));
    //cout << strings[i] << " " << strings[0] << " " << probs[i][0] << " " << endl;
    //cout << strings[0] << " " << strings[i] << " " << probs[0][i] << " " <<  endl;
    for (int j = 1; j < 3; ++j) {
      fst->AddArc(st, Arc(i, j, log(probs[i][j]), st));   
      //cout << strings[i] << " " << strings[j] << " " << probs[i][j] << " " <<  endl;
    }
  }
}


template<class Arc>
void BuildLinearChainFst(const string& str, VectorFst<Arc>* ofst, SymbolTable* syms) {
  typedef typename Arc::Label L;
  typedef typename Arc::StateId S;
  typedef typename Arc::Weight W;
  char* cstr = new char[str.size() + 1];
  strcpy(cstr, str.c_str());
  vector<char*> fields;
  SplitToVector(cstr, " ", &fields, true);
  S s = ofst->AddState();
  ofst->SetStart(s);
  for (size_t i = 0; i != fields.size(); ++i) {
    S d = ofst->AddState();
    L l = syms->AddSymbol(fields[i]);
    ofst->AddArc(s, Arc(l, l, W::One(), d));
    s = d;
  }
  ofst->SetFinal(s, W::One());
  delete[] cstr;
}



int mainfst() {
  LogVectorFst E;
  InitEditFst(&E, 3);
  SymbolTable syms("syms");
  syms.AddSymbol("<eps>", 0);

  LogVectorFst X;
  BuildLinearChainFst("<eps> H T", &X, &syms);
  LogVectorFst Y;
  BuildLinearChainFst("<eps> T T T", &Y, &syms);

  LogVectorFst XE;
  LogVectorFst XEY;

  Compose(X, E, &XE);
  Compose(XE, Y, &XEY);

  E.Write("E.ofst");
  XE.Write("XE.ofst");
  XEY.Write("XEY.ofst");
  
  typedef LogArc::Weight W;
  vector<W> alpha;
  vector<W> beta;
  ShortestDistance(XEY, &alpha);
  ShortestDistance(XEY, &beta, true);


  cout << "Forward probs" << endl;
  for (int i = 0; i < alpha.size(); ++i) 
    cout << i << " " << alpha[i].Value() << " " << exp(-alpha[i].Value()) << endl;

  cout << endl <<  "Backward probs" << endl;
  for (int i = 0; i < beta.size(); ++i) 
    cout << i << " " <<beta[i].Value() << " " << exp(-beta[i].Value()) << endl;
  return 0;
}


template<class Arc>
class EMTrainer {
  typedef typename Arc::Weight Weight;
  typedef typename Arc::Label Label;
  typedef typename Arc::StateId StateId;
  vector< vector<Weight> > alphas_;
  vector< vector<Weight> > betas_;
  vector< vector<Weight> > deltas_; //Probability distribution
  vector< vector<Weight> > deltas2_;
  Weight delta_final_;

  const vector< vector<int> >& iseqs_;
  const vector< vector<int> >& oseqs_;
  const SymbolTable& isyms_;
  const SymbolTable& osyms_;

 public:
  EMTrainer(const vector< vector<int> >& iseqs, 
            const vector< vector<int> >& oseqs,
            const SymbolTable& isyms, const SymbolTable& osyms):
    iseqs_(iseqs), oseqs_(oseqs), isyms_(isyms), osyms_(osyms) {

    float p = log(isyms.NumSymbols() * osyms.NumSymbols());
    LOG(INFO) << "p " << p;
    AllocMatrix(isyms.NumSymbols(), osyms.NumSymbols(), &deltas_, p);
    AllocMatrix(isyms.NumSymbols(), osyms.NumSymbols(), &deltas2_, Weight::Zero());
    delta_final_ = p;
    deltas_[0][0] = 0;

    size_t imax = 0;
    for (size_t i = 0; i != iseqs_.size(); ++i) {
      imax = max(imax, iseqs_[i].size());
      for (size_t j = 0; j < iseqs_[i].size(); ++j) 
        if (iseqs_[i][j] >= isyms_.NumSymbols())
          LOG(FATAL) << "Input symbol check failures " << i << " " << j << " " << iseqs_[i][j];
    }

    size_t omax = 0;
    for (size_t i = 0; i != oseqs_.size(); ++i) {
      omax = max(omax, oseqs_[i].size());      
      for (size_t j = 0; j < oseqs_[i].size(); ++j) 
        if (oseqs_[i][j] >= osyms_.NumSymbols())
          LOG(FATAL) << "Input symbol check failures " << i << " " << j << " " << oseqs_[i][j];
    }

    LOG(INFO) << "imax=" << imax << " omax=" << omax;
    AllocMatrix(imax, omax, &alphas_, Weight::Zero());
    AllocMatrix(imax, omax, &betas_, Weight::Zero());
    //LOG(INFO) << "Alphas";
    //PrintMatrix(alphas_);
    //LOG(INFO) << "Betas";
    //PrintMatrix(betas_);
    //LOG(INFO) << "Deltas";
    //PrintMatrix(deltas_);    
  }

  void AllocMatrix(size_t M, size_t N, vector<vector<Weight> >* matrix, Weight w) {
    matrix->clear();
    vector<Weight> row;
    for (size_t i = 0; i != N; ++i)
      row.push_back(w);
    for (size_t i = 0; i != M; ++i) 
      matrix->push_back(row);
  }

  void PrintMatrix(const vector< vector<Weight> >& matrix) {
    stringstream ss;
    for (int i = 0; i != matrix.size(); ++i) {
      const vector<Weight>& row = matrix[i];
      for (int j = 0; j != row.size(); ++j) {
        ss << row[j];
        if (j < row.size() - 1)
          ss << "\t";
        else
          ss << "\n";
      }
    }
    cerr << ss.str();
  }
  
  void PrintMatrixReal(const vector< vector<Weight> >& matrix) {
    stringstream ss;
    for (int i = 0; i != matrix.size(); ++i) {
      const vector<Weight>& row = matrix[i];
      for (int j = 0; j != row.size(); ++j) {
        ss << exp(-row[j].Value());
        if (j < row.size() - 1)
          ss << "\t";
        else
          ss << "\n";
      }
    }
    cerr << ss.str();
  }


  void ResetMatrix(vector< vector<Weight> >* matrix) {
    vector< vector<Weight> >& m = *matrix;
    for (size_t i = 0; i != m.size(); ++i) 
      for (size_t j = 0; j != m[i].size(); ++j)
        m[i][j] = Weight::Zero();
  }
 
  Weight Alphas(const vector<int>& is, const vector<int>& os, MutableFst<Arc>* fst = 0, SymbolTable* ssyms = 0) {
    ResetMatrix(&alphas_);
    alphas_[0][0] = Weight::One();
    size_t M = is.size();
    size_t N = os.size();
    for (size_t m = 0; m != M; ++m)
      for(size_t n = 0; n != N; ++n) {
        Weight& a =  alphas_[m][n];
        int x = is[m];
        int y = os[n];
        if (x >= deltas_.size())
          LOG(FATAL) << "X out of bounds " << x << " " << deltas_.size() << " " << m;
        if (y >= deltas_[x].size()) {
          for (int z = 0; z != N; ++z)
            cout << os[z] << " ";
          cout << endl;
          LOG(FATAL) << "Y out of bounds " << y << " " << deltas_[x].size() << " " << n;
        }
        if (m > 0) {
          Weight d = Times(alphas_[m - 1][n], deltas_[x][0]);
          a = Plus(a, d);
        }
        if (n > 0) {
          Weight i = Times(alphas_[m][n - 1], deltas_[0][y]);
          a = Plus(a, i);
        }
        if (m > 0 && n > 0) {
          Weight s = Times(alphas_[m - 1][n - 1], deltas_[x][y]);
          a = Plus(a, s);
        }
      }
    //PrintMatrix(alphas_);
    //Add the final accept cost
    Weight& final = alphas_[M - 1][N - 1];
    final = Times(final, delta_final_);

    if (fst) {
      fst->SetInputSymbols(&isyms_);
      fst->SetOutputSymbols(&osyms_);
      for (size_t i = 0; i != M * N; ++i) 
        fst->AddState();
      fst->SetStart(0);
      fst->SetFinal(M * N - 1, delta_final_);
      for (size_t m = 0; m != M; ++m) {
        for(size_t n = 0; n != N; ++n) {
          int d = m * N + n;
          int x = is[m];
          int y = os[n];
          if (m > 0) {
            int s = (m - 1) * N + n;
            fst->AddArc(s, LogArc(x, y, deltas_[x][0], d));
          }
          if (n > 0) {
            int s = m * N + n - 1;
            fst->AddArc(s, LogArc(x, y, deltas_[0][y], d));
          }
          if (m > 0 && n > 0) {
            int s = (m - 1) * N + n - 1;
            fst->AddArc(s, LogArc(x, y, deltas_[x][y], d));
          }
        }
      }
    }
    return final;
  }

  Weight Betas(const vector<int>& is, const vector<int>& os) {
    ResetMatrix(&betas_);
    size_t M = is.size() - 1;
    size_t N = os.size() - 1;
    betas_[M][N] = delta_final_;
    for (int m = M; m >= 0; --m)
      for(int n = N; n >= 0; --n) {
        Weight& b = betas_[m][n];
        int x = is[m + 1];
        int y = os[n + 1];
        if (m < M) {
          Weight d = Times(betas_[m + 1][n], deltas_[x][0]);
          b = Plus(b, d);
        }        
        if (n < N) {
          Weight i = Times(betas_[m][n + 1], deltas_[0][y]);
          b = Plus(b, i);
        }

        if (m < M  && n < N) {
          Weight s = Times(betas_[m + 1][n + 1], deltas_[x][y]);
          b = Plus(b, s);
        }
      }
    //PrintMatrix(betas_);
    return betas_[0][0];
  }

  Weight Gamma(Weight w1, Weight w2, Weight w3, Weight w4) {
    //LOG(INFO) << w1 << "," << w2 << "," << w3 << "," << w4;
    return Divide(Times(Times(w1, w2), w3), w4);
  }
  //Performs a round of EM training and returns the total log probability
  //over the test set
  double EMStep() {
    ResetMatrix(&deltas2_);
    Weight sum = Weight::One();
    double aveprob = 0.0f;
    for (size_t i = 0; i != iseqs_.size(); ++i) {
      //LOG(INFO) << i << endl;
      const vector<int>& is = iseqs_[i];
      const vector<int>& os = oseqs_[i];
      LogVectorFst fst;
      Weight fwd = Alphas(is, os, &fst);
      fst.Write("alpha.fst");
      Weight bwd = Betas(is, os);
      aveprob += fwd.Value();
      for (int m = 0; m != is.size(); ++m)
        for (int n = 0; n != os.size(); ++n) {
          int x = is[m];
          int y = os[n];
          if (x >= deltas_.size())
            LOG(FATAL) << "X out of bounds " << x << " " << deltas_.size() << " " << m;
          if (y >= deltas_[x].size())
            LOG(FATAL) << "Y out of bounds " << y << " " << deltas_[x].size() << " " << n;
          if (m > 0) {
            Weight w = Gamma(alphas_[m - 1][n], deltas_[x][0], betas_[m][n], betas_[0][0]);
            deltas2_[x][0] = Plus(deltas2_[x][0], w);
            sum = Plus(sum, w);
          }
          if (n > 0) {
            Weight w = Gamma(alphas_[m][n - 1], deltas_[0][y], betas_[m][n], betas_[0][0]);
            deltas2_[0][y] = Plus(deltas2_[0][y], w);
            sum = Plus(sum, w);
          }
          if (m > 0 && n > 0) {
            Weight w = Gamma(alphas_[m - 1][n - 1], deltas_[x][y], betas_[m][n], betas_[0][0]);
            deltas2_[x][y] = Plus(deltas2_[x][y], w);
            sum = Plus(sum, w);
          }
        }
    }
    for (int m = 0; m != isyms_.NumSymbols(); ++m)
      for (int n = 0; n != osyms_.NumSymbols(); ++n) {
        deltas2_[m][n] = Divide(deltas2_[m][n], sum);
      }
    Weight final = Weight::One();
    deltas_ = deltas2_;
    delta_final_ = Divide(final, sum);
    deltas_[0][0] = Weight::One();
    return aveprob / iseqs_.size();
  }

  void PrintDeltas() {
    PrintMatrixReal(deltas_);
    PrintMatrixReal(deltas_);
  }
  //Construct a new Fst based on the current table
  void GetFst(MutableFst<Arc>* fst) {
    fst->DeleteStates();
    fst->AddState();
  }

  void Align(ofstream& ofs) {
    vector< vector<float> > distances;
    vector< vector<int> > back;
    int omax = alphas_[0].size();
    vector<float> drow;
    drow.resize(omax, 1000000);
    vector<int> brow;
    brow.resize(omax, -1);

    for (int i = 0; i != alphas_.size(); ++i) {
      distances.push_back(drow);
      back.push_back(brow);
    }

    for (int i = 0; i != iseqs_.size(); ++i ) {
      const vector<int>& is = iseqs_[i];
      const vector<int>& os = oseqs_[i];
      size_t M = is.size();
      size_t N = os.size();
      for (size_t m = 0; m != M; ++m)
        for(size_t n = 0; n != N; ++n)
          distances[m][n] = 100000;
      
      distances[0][0] = 0;
      for (size_t m = 0; m != M; ++m)
        for(size_t n = 0; n != N; ++n) {
          int x = is[m];
          int y = os[n];
          float& f = distances[m][n];
          int& b = back[m][n];
          if (m > 0) {
            float d =  distances[m - 1][n] + deltas_[x][0].Value();
            if (d < f) {
              f = d;
              b = 0;
            }
          }
          if (n > 0) {
            float i = distances[m][n - 1] + deltas_[0][y].Value();
            if (i < f) {
              f = i;
              b = 1;
            }
          }
          if (m > 0 && n > 0) {
            float s = distances[m - 1][n - 1] + deltas_[x][y].Value();
            //LOG(INFO) << s << " " << f;
            //cin.get();
            if (s < f) {
              f = s;
              b = 2;
            }
          }
       }
      
      /*
      for (size_t m = 0; m != M; ++m)
        for(size_t n = 0; n != N; ++n) {
          cout << m << " "  << n << " " << distances[m][n] << " " << back[m][n] << endl;
        }
      */
       size_t m = M - 1;
       size_t n = N - 1;
       //cout << "best alignment " << distances[m][n] << endl;
       int b = back[m][n];
       vector< pair<string,string> > alignment;
       const char* space = "\xe3\x80\x80";
       while (b != -1) {
         //cout << b << " : ";
         string s = space;
         string t = space;
         switch(b) {
           case 2: { s = isyms_.Find(is[m]); t = osyms_.Find(os[n]); --m; --n; b = back[m][n];  break; }
           case 1: { t = osyms_.Find(os[n]); --n; b = back[m][n];  break; }
           case 0: { s = isyms_.Find(is[m]); --m; b = back[m][n];  break; }
         }
         alignment.push_back(pair<string,string>(s,t));
       }
       //cout << endl;
       stringstream sss;
       stringstream sst;
       for (int i = alignment.size() - 1; i >= 0; --i) {
         pair<string,string> p = alignment[i];
         sss << p.first;
         sst << p.second;
       }
       ofs << sss.str() << endl;
       ofs << sst.str() << endl << endl;
       //backtrack
       //break;
    }     
  }
};

void SplitStringToVector(const std::string &full, const char *delim,
                         bool omit_empty_strings,
                         std::vector<std::string> *out) {
  size_t start = 0, found = 0, end = full.size();
  out->clear();
  while (found != std::string::npos) {
    found = full.find_first_of(delim, start);
    // start != end condition is for when the delimiter is at the end
    if (!omit_empty_strings || (found != start && start != end))
      out->push_back(full.substr(start, found - start));
    start = found + 1;
  }
}


DEFINE_int32(num_iter, 1, "Number of training iterations");

int main(int argc, char** argv) {
  
  string usage = "Train a stochastic edit distance Fst.\n\n  Usage: ";
  usage += argv[0];
  usage += " [in.train.text [out.align.text]]\n";

  std::set_new_handler(FailedNewHandler);
  SET_FLAGS(usage.c_str(), &argc, &argv, true);

  if (argc > 3) {
    ShowUsage();
    return 1;
  }

  SymbolTable isyms("isyms");
  SymbolTable osyms("osyms");
  isyms.AddSymbol("<eps>", 0);
  osyms.AddSymbol("<eps>", 0);
  
  vector< vector<int> > iseqs;
  vector< vector<int> > oseqs;

  string ip = argc >= 2 ? argv[1] : "/dev/stdin";
  LOG(INFO) << "Input file : " << ip;
  ifstream ifs(ip.c_str());
  if (! ifs.is_open())
    LOG(FATAL) << "Failed to open training data : ";

  string op = argc >= 2 ? argv[2] : "/dev/stdout";
  LOG(INFO) << "Output file : " << op;
  ofstream ofs(op.c_str());
  if (! ofs.is_open())
    LOG(FATAL) << "Failed to open alignment file : ";


  //Read the corpus line by line and
  //add an epsilon transition to the start of each i/o seq
  //this corresponds to an insertion or deletion
  for (string s; getline(ifs, s); ) {
    vector<string> ios;
    vector<string> is;
    vector<string> os;
    SplitStringToVector(s, "\t", true, &ios);
    if (ios.size() != 2)
      LOG(FATAL) << "Incorrect number of fields " << s;
    SplitStringToVector(ios[0], " ", true, &is);
    SplitStringToVector(ios[1], " ", true, &os);

    if (is.empty())
      continue;

    if (os.empty())
      continue;

    vector<int> i;
    i.push_back(0);
    for (size_t n = 0; n != is.size(); ++n)
      i.push_back(isyms.AddSymbol(is[n]));
    iseqs.push_back(i);
    
    vector<int> o;
    o.push_back(0);
    for (size_t n = 0; n != os.size(); ++n)
      o.push_back(osyms.AddSymbol(os[n]));
    oseqs.push_back(o);
  }

  LOG(INFO) << "Read " << iseqs.size() << " training pairs";
  LOG(INFO) << "Num input symbols " << isyms.NumSymbols();
  LOG(INFO) << "Num output symbols " << osyms.NumSymbols();
  //Perform EM training for n iterations
  EMTrainer<LogArc> trainer(iseqs, oseqs, isyms, osyms);
  LOG(INFO) << "Beginning training procedure";
  for (int n = 0; n != FLAGS_num_iter; ++n) {
    double d = trainer.EMStep();
    LOG(INFO) << "Done " << n << " " << d;
    //trainer.PrintDeltas();
  }
  LOG(INFO) << "Aligning training data";
  trainer.Align(ofs);
  LOG(INFO) << "Finished everything";
  //Write the trained Fst to disk
  return 0;
}
