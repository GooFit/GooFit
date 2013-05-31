#include "Variable.hh"
#include <cmath> 

//std::map<std::string, Variable*> Variable::nameToVariableMap;

Variable::Variable (std::string n) 
  : name(n) 
  , numbins(100)
  , index(-1) 
  , fixed(false)
  , blind(0)
{
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
  , blind(0)
{
}

Variable::Variable (std::string n, fptype dn, fptype up) 
  : name(n)
  , upperlimit(up)
  , lowerlimit(dn)
  , numbins(100)
  , index(-1)
  , fixed(false)
  , blind(0)
{
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
  , blind(0)
{
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
  , blind(0)
{
}

Variable::~Variable () {
}

