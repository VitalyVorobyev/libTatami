/** Copyright 2016 Vitaly Vorobyev
 ** @file icpvevent.cpp
 **
 ** @brief This message displayed in Doxygen Files index
 **
 ** @author Vitaly Vorobyev
 ** Contact: vit.vorobiev@gmail.com
 **
 **/

#include <fstream>
#include <iostream>
#include <sstream>
#include <exception>

#include "icpvevent.h"

namespace libTatami {

using std::istringstream;
using std::cout;
using std::cerr;
using std::endl;
using std::getline;
using std::ifstream;
using std::vector;
using std::string;

ICPVEvt::ICPVEvt(const string& fname) {
    ReadStructure(fname);
}

ICPVEvt::ICPVEvt(const ICPVEvt& ev) {
    Set(ev.IVars(), ev.DVars());
}

ICPVEvt::ICPVEvt(const ivarvec& IV, const dvarvec& DV) {
    Set(IV, DV);
}

ICPVEvt& ICPVEvt::operator=(const ICPVEvt& vt) {
    Set(vt.IVars(), vt.DVars());
    return *this;
}

int ICPVEvt::ReadStructure(const string& fname) {
    ifstream ifile(fname.c_str(), ifstream::in);
    if (!ifile.is_open()) {
        cout << "ICPVEvt::ReadStructure: can't open file " << fname << endl;
        return -1;
    } else {
        cout << "Getting event structure from file " << fname << endl;
    }
    m_IVars.clear(); m_DVars.clear();
    string line;
    while (getline(ifile, line)) {
        string name, type, dscr;
        istringstream out; out.str(line);
        out >> name; out >> type; out >> dscr;
        if      (type == string("D")) m_DVars.push_back(dvar(0, name));
        else if (type == string("I")) m_IVars.push_back(ivar(0, name));
        else
            cout << "ICPVEvt::ReadStructure: wrong var type " << type << endl;
    }
    PrintStructure();
    return 0;
}

void ICPVEvt::PrintStructure(void) const {
    cout << "Integer variables (" << m_IVars.size() << "):" << endl;
    for (unsigned i = 0; i < m_IVars.size(); i++)
        cout << " " << m_IVars[i].Name();
    cout << endl;
    cout << "Real variables (" << m_DVars.size() << "):" << endl;
    for (unsigned i = 0; i < m_DVars.size(); i++)
        cout << " " << m_DVars[i].Name();
    cout << endl;
    return;
}

string ICPVEvt::IName(uint32_t i) const {
    if (i < 0 || i >= m_IVars.size()) {
        cout << "ICPVEvt::IName: wrong var index " << i << endl;
        return 0;
    }
    return m_IVars[i].Name();
}

int ICPVEvt::IVar(uint32_t i) const {
    if (i < 0 || i >= m_IVars.size()) {
        cout << "ICPVEvt::IVar: wrong var index " << i << endl;
        return 0;
    }
    return m_IVars[i].Val();
}

string ICPVEvt::DName(uint32_t i) const {
    if (i < 0 || i >= m_DVars.size()) {
        cout << "ICPVEvt::DName: wrong var index " << i << endl;
        return 0;
    }
    return m_DVars[i].Name();
}

double ICPVEvt::DVar(uint32_t i) const {
    if (i < 0 || i >= m_DVars.size()) {
        cout << "ICPVEvt::DVar: wrong var index " << i << endl;
        return 0;
    }
    return m_DVars[i].Val();
}

int ICPVEvt::IVar(const string& name) const {
    return IVar(FindIVar(name));
}

double ICPVEvt::DVar(const string& name) const {
    return DVar(FindDVar(name));
}

void ICPVEvt::Set(const vector<ivar>& IV, const vector<dvar>& DV) {
    m_IVars = IV;
    m_DVars = DV;
}

void ICPVEvt::SetIVar(uint32_t i, uint32_t val) {
    if (i < 0 || i >= m_IVars.size()) {
        cout << "ICPVEvt::SetIVar: wrong var index " << i << endl;
        return;
    }
    m_IVars[i].SetVal(val);
}

void ICPVEvt::SetIVar(const string& name, uint32_t val) {
    m_IVars[FindIVar(name)].SetVal(val);
}

void ICPVEvt::SetDVar(uint32_t i, double val) {
    if (i < 0 || i >= m_DVars.size()) {
        cerr << "ICPVEvt::SetDVar: wrong var index " << i << endl;
        throw new std::runtime_error("");
    }
    m_DVars[i].SetVal(val);
}

void ICPVEvt::SetDVar(const string& name, double val) {
    m_DVars[FindDVar(name)].SetVal(val);
}

int ICPVEvt::FindIVar(const string& name) const {
    auto index = ivar::FindIndex(m_IVars, name);
    if (index == -1) {
        cout << "ICPVEvt::FindIVar: no variable with name " << name << endl;
        return 0;
    }
    return index;
}

int ICPVEvt::FindDVar(const string& name) const {
    auto index = dvar::FindIndex(m_DVars, name);
    if (index == -1) {
        cout << "ICPVEvt::FindDVar: no variable with name " << name << endl;
        return 0;
    }
    return index;
}

}  // namespace libTatami
