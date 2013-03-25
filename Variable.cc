  // Removed nameToVariableMap and related methods to make the class thread safe.
  // K. Tomko, May 11, 2012
#include "Variable.hh"
#include <cmath> 

//std::map<std::string, Variable*> Variable::nameToVariableMap;

Variable::Variable (std::string n) 
  : name(n) 
  , numbins(100)
  , index(-1) 
  , fixed(false)
{
  //  checkName(); 
} 

Variable::Variable (std::string n, fptype v) 
  : name(n)
  , value(v)
  , error(0.002) 
  , lowerlimit(v - 0.01)
  , upperlimit(v + 0.01)
  , numbins(100)
  , index(-1)
  , fixed(true)
{
  // checkName(); 
}

Variable::Variable (std::string n, fptype dn, fptype up) 
  : name(n)
  , upperlimit(up)
  , lowerlimit(dn)
  , numbins(100)
  , index(-1)
  , fixed(false)
{
  //  checkName(); 
}

Variable::Variable (std::string n, fptype v, fptype dn, fptype up) 
  : name(n)
  , value(v) 
  , error(0.1*(up-dn))
  , upperlimit(up)
  , lowerlimit(dn)
  , numbins(100)
  , index(-1)
  , fixed(false)
{
  //  checkName(); 
}

Variable::Variable (std::string n, fptype v, fptype e, fptype dn, fptype up) 
  : name(n)
  , value(v) 
  , error(e)
  , upperlimit(up)
  , lowerlimit(dn)
  , numbins(100)
  , index(-1)
  , fixed(false)
{
  //  checkName(); 
}

Variable::~Variable () {
  //  nameToVariableMap[name] = 0; 
}

//void Variable::checkName () {
//  if (nameToVariableMap[name]) {
//    std::cout << "Error: Variable with name " 
//	      << name
//	      << " already exists, aborting.\n";
//    assert(false);
//  }
//  nameToVariableMap[name] = this; 
//}

//Variable* Variable::getVariable (std::string n) {
  /*if (!nameToVariableMap[n]) {
    std::cout << "Error: No variable named "
	      << n
	      << ", aborting.\n";
    assert(false);
    }*/
//  return nameToVariableMap[n]; 
//}
