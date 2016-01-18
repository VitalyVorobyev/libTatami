#include "ResConst.h"

#include <fstream>

using namespace std;

ResConst::ResConst(const std::string& fname){
  LoadParameters(fname);
  return;
}

int ResConst::LoadParameters(const std::string& fname){
  ifstream ifile(fname.c_str(),ofstream::in);
  if(!ifile.is_open()){
    cout << "ResConst::LoadParameters: can't open file " << fname << endl;
    return -1;
  } else{
    cout << "Getting dt resolution from file " << fname << endl;
  }
  string line,name;
  double val,errp,errn;
  char namech[15];
  int counter = 0;
  for(int i=0; i<13; i++){
    getline(ifile,line);
    if(4 != sscanf(line.c_str(),"%s = %lf %lf %lf",namech,&val,&errp,&errn)){
      if(2 != sscanf(line.c_str(),"%s = %lf",namech,&val)){
        cout << "ResConst::LoadParameters: can't parse string " << line << endl;
        continue;
      }
    }
    name = string(namech);
    cout << name << " " << val << endl;
    if(name == string("Srec0"))             { m_Srec[0]           = val; counter++; continue;}
    if(name == string("Srec1"))             { m_Srec[1]           = val; counter++; continue;}
    if(name == string("Sasc0"))             { m_Sasc[0]           = val; counter++; continue;}
    if(name == string("Sasc1"))             { m_Sasc[1]           = val; counter++; continue;}
    if(name == string("Snp_global"))        { m_Snp_global        = val; counter++; continue;}
    if(name == string("Smn_rec"))           { m_Smn_rec           = val; counter++; continue;}
    if(name == string("Stl_rec"))           { m_Stl_rec           = val; counter++; continue;}
    if(name == string("ftl_rec"))           { m_ftl_rec           = val; counter++; continue;}
    if(name == string("Smn_asc"))           { m_Smn_asc           = val; counter++; continue;}
    if(name == string("Stl_asc"))           { m_Stl_asc           = val; counter++; continue;}
    if(name == string("ftl_asc"))           { m_ftl_asc           = val; counter++; continue;}

    if(name == string("fd_np_sgl"))         { m_fd_np_sgl         = val; counter++; continue;}
    if(name == string("tau_np_p_sgl0"))     { m_tau_np_p_sgl[0]   = val; counter++; continue;}
    if(name == string("tau_np_p_sgl1"))     { m_tau_np_p_sgl[1]   = val; counter++; continue;}
    if(name == string("tau_np_n_sgl0"))     { m_tau_np_n_sgl[0]   = val; counter++; continue;}
    if(name == string("tau_np_n_sgl1"))     { m_tau_np_n_sgl[1]   = val; counter++; continue;}
    if(name == string("fd_np_mlt0"))        { m_fd_np_mlt[0]      = val; counter++; continue;}
    if(name == string("fd_np_mlt1"))        { m_fd_np_mlt[1]      = val; counter++; continue;}
    if(name == string("fd_np_st_mlt"))      { m_fd_np_st_mlt      = val; counter++; continue;}
    if(name == string("fd_np_xi_mlt"))      { m_fd_np_xi_mlt      = val; counter++; continue;}
    if(name == string("fd_np_stxi_mlt"))    { m_fd_np_stxi_mlt    = val; counter++; continue;}
    if(name == string("fn_np_mlt"))         { m_fn_np_mlt         = val; counter++; continue;}
    if(name == string("tau_np_p_mlt0"))     { m_tau_np_p_mlt[0]   = val; counter++; continue;}
    if(name == string("tau_np_p_mlt1"))     { m_tau_np_p_mlt[1]   = val; counter++; continue;}
    if(name == string("tau_np_p_xi_mlt"))   { m_tau_np_p_xi_mlt   = val; counter++; continue;}
    if(name == string("tau_np_p_stxi_mlt")) { m_tau_np_p_stxi_mlt = val; counter++; continue;}
    if(name == string("tau_np_n_mlt0"))     { m_tau_np_n_mlt[0]   = val; counter++; continue;}
    if(name == string("tau_np_n_mlt1"))     { m_tau_np_n_mlt[1]   = val; counter++; continue;}
    if(name == string("tau_np_n_xi_mlt"))   { m_tau_np_n_xi_mlt   = val; counter++; continue;}
    if(name == string("tau_np_n_stxi_mlt")) { m_tau_np_n_stxi_mlt = val; counter++; continue;}

    if(name == string("sigma_ol"))          { m_sigma_ol          = val; counter++; continue;}
    if(name == string("f_ol_sgl"))          { m_f_ol_sgl          = val; counter++; continue;}
    if(name == string("f_ol_mul"))          { m_f_ol_mul          = val; counter++; continue;}

    if(name == string("rnp_kink_xi"))       { m_rnp_kink_xi       = val; counter++; continue;}
    if(name == string("rnp_kink_st"))       { m_rnp_kink_st       = val; counter++; continue;}
  }
  return counter;
}

