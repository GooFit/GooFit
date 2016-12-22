#include "goofit/DataSet.h" 
#include <cstdlib> 
#include <climits> 

DataSet::DataSet (Variable* var, string n)  
  : numEventsAdded(0)
  , name(n)
{
  variables.push_back(var); 
  if (n == "") generateName();
}

DataSet::DataSet (std::vector<Variable*>& vars, string n)  
  : numEventsAdded(0)
  , name(n)
{
  for (std::vector<Variable*>::iterator v = vars.begin(); v != vars.end(); ++v) {
    variables.push_back(*v); 
  }
  if (n == "") generateName();
}

DataSet::DataSet (std::set<Variable*>& vars, string n)  
  : numEventsAdded(0)
  , name(n)
{
  variables.resize(vars.size()); 
  
  for (std::set<Variable*>::iterator v = vars.begin(); v != vars.end(); ++v) {
    variables[(*v)->index] = (*v); 
  }
  if (n == "") generateName();
}

DataSet::~DataSet () {
  variables.clear(); 
}

void DataSet::addEvent () {
  vector<fptype> vals = getCurrentValues(); 
  addEventVector(vals); 
}

void DataSet::addWeightedEvent (fptype weight) {
  vector<fptype> vals = getCurrentValues(); 
  addEventVector(vals, weight); 
}

void DataSet::addEvent (fptype val) {
  // Helper method to avoid the user having to wrap 
  // every one-dimensional event in a vector. 
  assert(1 == variables.size()); 

  std::vector<fptype> helper;
  helper.push_back(val); 
  addEventVector(helper); 
}

vector<fptype> DataSet::getCurrentValues () const {
  vector<fptype> values;
  for (varConstIt v = varsBegin(); v != varsEnd(); ++v) {
    values.push_back((*v)->value);
  }
  return values; 
}

void DataSet::getVariables (std::vector<Variable*>& vars) {
  for (std::vector<Variable*>::iterator v = variables.begin(); v != variables.end(); ++v) {
    vars.push_back(*v);
  }
}

unsigned int DataSet::indexOfVariable (Variable* var) const {
  for (unsigned int i = 0; i < variables.size(); ++i) {
    if (var != variables[i]) continue;
    return i; 
  }

  std::cout << "Error: Attempted to get index of variable "
	    << var->name
	    << " in DataSet of ";
  for (varConstIt v = varsBegin(); v != varsEnd(); ++v) {
    std::cout << "\n  " << (*v)->name << std::endl;
  }
  std::cout << "\nAborting." << std::endl;
  assert(false); 
  abort(); 
  return 0; 
}

void DataSet::generateName () {
  // Create default name as list of variables.
  if (name != "") return; 
  for (varConstIt v = varsBegin(); v != varsEnd(); ++v) {
    name += (*v)->name;
    varConstIt next = v; next++; if (next == varsEnd()) continue;
    name += ", "; 
  }
}
