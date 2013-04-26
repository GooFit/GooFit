#ifndef BINNED_DATASET_HH
#define BINNED_DATASET_HH

#include "DataSet.hh" 

class BinnedDataSet : public DataSet {
  // Class for rectangularly binned datasets - every bin the same size. 

public: 
  BinnedDataSet (Variable* var); 
  BinnedDataSet (std::vector<Variable*>& vars); 
  BinnedDataSet (std::set<Variable*>& vars); 
  ~BinnedDataSet (); 

  virtual void addEventVector (std::vector<fptype>& vals, fptype weight = 1);

  fptype getBinContent (unsigned int bin) const {return binvalues[bin];} 
  fptype getBinCenter (Variable* var, unsigned int bin) const;
  unsigned int getBinNumber () const; 
  fptype getBinVolume (unsigned int bin) const; 
  fptype getBinError  (unsigned int bin) const; 
  unsigned int getNumBins () const; 
  fptype getNumEvents () const; 

  void setBinContent (unsigned int bin, fptype value) {binvalues[bin] = value;}
  void setBinError (unsigned int bin, fptype error); 

private:
  void cacheNumBins (); 
  vector<unsigned int> convertValuesToBins (const vector<fptype>& vals) const; 
  unsigned int localToGlobal (std::vector<unsigned int>& locals) const; 
  void globalToLocal (std::vector<unsigned int>& locals, unsigned int global) const; 

  std::vector<fptype> binvalues; 
  std::vector<fptype> binerrors; 
  std::map<Variable*, int> cachedNumBins; // Store these numbers in case they change on the user end - vast confusion possible. 
};


#endif
