
#include "GeneralSkimmer.h"
#include "iostream"
#include "vector"

#include <getopt.h>
#include <climits>

#include <TChain.h>
#include <TString.h>



void processTChain( std::vector<TString> filenames, TString option, long n_events ); 

int main(int argc, char * argv[])
{
    // default values
    std::vector<TString> filenames;
    long n_events = LONG_MAX ;
    TString option = "";

    int c;
    while(1) {
      static struct option long_options[] = {
        {"input", required_argument, 0, 'i'},
        {"n_events", required_argument, 0, 'e'},
        {"output", required_argument, 0, 'o'}
      };

      int option_index = 0;
      c = getopt_long (argc, argv, "i:e:o:h",long_options, &option_index); 

      // Detect the end of the options
      if (c == -1)  break;

      switch(c) {
        case 'i':
          filenames.push_back(TString(optarg));
          break;
        case 'e':
          n_events = atoi(optarg);
          break;
        case 'o':
          option += "ofile="+TString(optarg) + ";"; 
          break;
        case 'h': // Intentional fall-through
        case '?':
          return 0;  
          break;
        default:
          abort();  
      }
    }


    processTChain( filenames, option, n_events);
    return 0;
}


void processTChain( std::vector<TString> filenames, TString option,  long n_events ) {

  
  // create TChain and add all filenames 
  TChain * tchain = new TChain("Tree");
  for ( std::size_t i = 0; i < filenames.size(); i++ ) {
    std::cout << " Adding file: " << filenames[i];
    // -1 option so file is checked (return 1 if file is correct or 0 otherwise)
    if (tchain->Add( filenames[i], -1 ) == 0) {
      std::cout << "\033[1;31m not found \033[0m\n" << std::endl; 
    } else {
      std::cout << "\033[1;32m found \033[0m\n" << std::endl; 
    } 
  }
  // get total number of entries in the chain
  long int nentries = tchain->GetEntries();
  std::cout << " Total number of entries in TChain is: " << nentries << std::endl;

  // create selector instance
  GeneralSkimmer * selector = new GeneralSkimmer();

  // process TChain
  tchain->Process(selector, option , n_events);

} 




