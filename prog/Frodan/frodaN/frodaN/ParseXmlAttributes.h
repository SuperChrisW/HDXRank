/*
 * ParseXmlAttributes.h
 *
 *  Created on: Jun 15, 2010
 *      Author: dwfarrel
 */

#ifndef PARSEXMLATTRIBUTES_H_
#define PARSEXMLATTRIBUTES_H_

#define TIXML_USE_TICPP
#include "ticpp.h"
#include <string>
#include <map>

class ParseXmlAttributes {
public:
  ParseXmlAttributes() {}
  virtual ~ParseXmlAttributes() {}
  void add( std::string name, std::string* destination ) {
    string_opts[name] = destination;
  }
  void add( std::string name, double* destination ) {
    double_opts[name] = destination;
  }
  void add( std::string name, int* destination ) {
    int_opts[name] = destination;
  }
  void add( std::string name, bool* destination ) {
    bool_opts[name] = destination;
  }
  void parse( ticpp::Element* e );
  void clear() {
    string_opts.clear();
    double_opts.clear();
    int_opts.clear();
    bool_opts.clear();
  }
private:
  std::map<std::string,std::string*> string_opts;
  std::map<std::string,double*> double_opts;
  std::map<std::string,int*> int_opts;
  std::map<std::string,bool*> bool_opts;
};

#endif /* PARSEXMLATTRIBUTES_H_ */
