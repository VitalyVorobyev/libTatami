#include "ResVar.h"

const bool RdetVar::RecSide = true;
const bool RdetVar::AscSide = false;

RdetVar::RdetVar(void):
 m_ntrk(1), m_sz(1.), m_chisq(1.), ndf_z(0)
{
}

RdetVar::RdetVar(const RdetVar& var){
  *this = var;
}

RdetVar& RdetVar::operator=(const RdetVar& var){
  m_ntrk  = var.m_ntrk;
  m_sz    = var.m_sz;
  m_chisq = var.m_chisq;
  ndf_z   = var.ndf_z;
  return *this;
}

int RdetVar::ReadVars(const ICPVEvt &evt, const bool type){
  if(type == RecSide){
    dz_rec      = evt.FindDVar("dz_rec");
    ntrk_rec    = evt.FindIVar("ntrk_rec");
    sz_rec      = evt.FindDVar("sz_rec");
    chisq_z_rec = evt.FindDVar("chisq_z_rec");
    ndf_z_rec   = evt.FindIVar("ndf_z_rec");
  } else{
    dz_rec      = evt.FindDVar("dz_asc");
    ntrk_rec    = evt.FindIVar("ntrk_asc");
    sz_rec      = evt.FindDVar("sz_asc");
    chisq_z_rec = evt.FindDVar("chisq_z_asc");
    ndf_z_rec   = evt.FindIVar("ndf_z_asc");
  }
  return 0;
}

