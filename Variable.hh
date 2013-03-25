#ifndef VARIABLE_HH
#define VARIABLE_HH

  // Removed nameToVariableMap and related methods to make the class thread safe.
  // K. Tomko, May 11, 2012

#include <string> 
#include <map> 
#include <iostream> 
#include <cassert> 
#include "GlobalCudaDefines.hh"

struct Variable { 
  Variable (std::string n); 
  Variable (std::string n, fptype val); 
  Variable (std::string n, fptype dn, fptype up);
  Variable (std::string n, fptype v, fptype dn, fptype up);
  Variable (std::string n, fptype v, fptype e, fptype dn, fptype up);
  ~Variable (); 

  // static Variable* getVariable (std::string n); 
  int getIndex () const {return index;}

  std::string name;
  fptype value;
  fptype error;
  fptype upperlimit;
  fptype lowerlimit;
  int numbins; 
  int index; 
  bool fixed; 

  //private:
  //static std::map<std::string, Variable*> nameToVariableMap; 
  //void checkName ();
}; 

#endif
