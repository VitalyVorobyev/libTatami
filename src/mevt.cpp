/** Copyright 2016 Vitaly Vorobyev
 ** @file mevt.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include "mevt.h"

#include <iostream>
#include <fstream>
#include <sstream>

using std::cout;
using std::endl;
using std::ifstream;
using std::istringstream;
using std::getline;
using std::string;

namespace libTatami {

MEvt::MEvt(const string& fname) {
    ReadStructure(fname);
}

int MEvt::ReadStructure(const string& fname) {
    ifstream ifile(fname.c_str(), ifstream::in);
    if (!ifile.is_open()) {
        cout << "MEvt::ReadStructure: can't open file " << fname << endl;
        return -1;
    } else {
        cout << "Getting event structure from file " << fname << endl;
    }
    m_ivars.clear(); m_dvars.clear();
    string line, name, type, dscr;
    istringstream out;
    while (getline(ifile, line)) {
        out.str(line); out >> name; out >> type; out >> dscr;
        if (type == string("D")) m_dvars.insert(ivar(name, 0));
        else if (type == string("I")) m_ivars.insert(ivar(name, 0));
        else
            cout << "MEvt::ReadStructure: wrong var type " << type << endl;
    }
    PrintStructure();
    return 0;
}

void MEvt::PrintStructure(void) const {
    cout << "Integer variables (" << m_ivars.size() << "):" << endl;
    for (auto const& x : m_ivars) cout << " " << x.first;
    cout << endl;
    cout << "Real variables (" << m_dvars.size() << "):" << endl;
    for (auto const& x : m_dvars) cout << " " << x.first;
    cout << endl;
}

}  // namespace libTatami
