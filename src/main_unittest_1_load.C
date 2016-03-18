#define CATCH_CONFIG_MAIN

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cstdio>


#include "size.h"
#include "geauxdock.h"
#include "util.h"
#include "load.h"

#include <catch/catch.hpp>

TEST_CASE( "Load", "[ligand]" )
{
  Ligand0 *lig0 = new Ligand0[300];

  LigandFile *lig_file = new LigandFile;
  lig_file->molid = "MOLID";
  lig_file->path = "../data/test/ZINC00062388.sdf";

  loadLigand (lig_file, lig0);

  delete[]lig0;
}

/*
TEST_CASE( "LoadLHM", "[force field input]" )
{

  LhmFile *lhm_file = new LhmFile;
  Psp0 *psp = new Psp0;
  Kde0 *kde = new Kde0;
  Mcs0 *mcs = new Mcs0;

  lhm_file->path = "../data/test/1b9vA_3.ff";
  lhm_file->ligand_id = "ZINC49076157";


  loadLHM(lhm_file, psp, kde, mcs);

  delete mcs;
  delete psp;
  delete kde;
  delete lhm_file;
}
*/
