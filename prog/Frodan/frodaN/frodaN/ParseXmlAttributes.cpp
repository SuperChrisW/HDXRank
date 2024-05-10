/*
 * ParseXmlAttributes.cpp
 *
 *  Created on: Jun 15, 2010
 *      Author: dwfarrel
 */

#include "ParseXmlAttributes.h"
#include <sstream>
#include <iostream>

using namespace std;

void ParseXmlAttributes::parse( ticpp::Element* e ) {
  ticpp::Iterator<ticpp::Attribute> a;
  string aname;
  for ( a = a.begin( e ); a != a.end(); a++ ) {
    aname = a->Name();

    map<string,string*>::iterator s = string_opts.find( aname );
    if ( s != string_opts.end() ) {
      a->GetValue( s->second );
      continue;
    }

    map<string,double*>::iterator d = double_opts.find( aname );
    if ( d != double_opts.end() ) {
      a->GetValue( d->second );
      continue;
    }

    map<string,int*>::iterator i = int_opts.find( aname );
    if ( i != int_opts.end() ) {
      a->GetValue( i->second );
      continue;
    }

    map<string,bool*>::iterator b = bool_opts.find( aname );
    if ( b != bool_opts.end() ) {
      string temp;
      a->GetValue( &temp );

      if ( temp == "1" || temp == "true" || temp == "on" ) *b->second = true;
      else if ( temp == "0" || temp == "false" || temp == "off" ) *b->second = false;
      else {
        ostringstream o;
        o << "Attribute " << a->Name() << " not boolean, line " << a->Row() << '\n';
        o << "Allowed values are 0, 1, true, false, on, off (lowercase only)" << '\n';
        throw ticpp::Exception( o.str() );
      }
      continue;
    }

    ostringstream o;
    o << "Attribute " << a->Name() << " not recognized, line " << a->Row() << endl;
    throw ticpp::Exception( o.str() );

  }
}
