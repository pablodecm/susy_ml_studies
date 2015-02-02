
#include "GeneralSkimmer.h"
#include "iostream"

#include <TChain.h>

int main(int argc, char * argv[])
{
    std::cout << "The GeneralSkimmer will be applied to one ntuple" << std::endl;
    GeneralSkimmer * selector = new GeneralSkimmer();
    TChain * tchain = new TChain("Tree");
    tchain->Add("/gpfs/csic_projects/tier3data/MC_Summer12_53X/SUSY/SMS/Tree_T2tt_150to350LSP0to250_0.root");
    tchain->Process(selector, "", 50000); //only process 50000 entries
    return 0;
}
