#ifndef DATASET_HH
#define DATASET_HH

#include "Variable.hh"
#include "GlobalCudaDefines.hh"
#include <vector> 
#include <set> 

typedef std::vector<Variable*>::const_iterator varConstIt; 
typedef std::vector<Variable*>::const_reverse_iterator varConstRIt; 

class DataSet {
public:
  DataSet (Variable* var); 
  DataSet (std::vector<Variable*>& vars); 
  DataSet (std::set<Variable*>& vars); 
  ~DataSet (); 

  void addEvent ();
  void addEvent (fptype val);
  virtual void addEventVector (std::vector<fptype>& vals, fptype weight = 1) = 0; 
  void addWeightedEvent (fptype weight); 

  varConstIt varsBegin () const {return variables.begin();}
  varConstIt varsEnd () const {return variables.end();}
  void getVariables (std::vector<Variable*>& vars); 

  varConstRIt varsRBegin () const {return variables.rbegin();}
  varConstRIt varsREnd () const {return variables.rend();}
  int numVariables () const {return variables.size();} 
  int numEvents () const {return numEventsAdded;}

protected:
  vector<fptype> getCurrentValues () const; 
  unsigned int indexOfVariable (Variable* var) const;
  int numEventsAdded; 

private: 
  std::vector<Variable*> variables; 
}; 



#endif
