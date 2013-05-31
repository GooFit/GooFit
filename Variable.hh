#ifndef VARIABLE_HH
#define VARIABLE_HH

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

  int getIndex () const {return index;}

  std::string name;
  fptype value;
  fptype error;
  fptype upperlimit;
  fptype lowerlimit;
  int numbins; 
  int index; 
  bool fixed; 
  fptype blind; 
}; 

#endif
