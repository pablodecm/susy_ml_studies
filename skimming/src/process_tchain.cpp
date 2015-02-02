
#include "GeneralSkimmer.h"
#include "iostream"
#include "vector"

#include <getopt.h>
#include <climits>

#include <TChain.h>
#include <TString.h>



void processTChain( std::vector<TString> filenames, long n_events ); 

int main(int argc, char * argv[])
{
    // default values
    std::vector<TString> filenames;
    long n_events = LONG_MAX ;

    int c;
    while(1) {
      static struct option long_options[] = {
        {"input", required_argument, 0, 'i'},
        {"n_events", required_argument, 0, 'e'}
      };

      int option_index = 0;
      c = getopt_long (argc, argv, "i:e:h",long_options, &option_index); 

      // Detect the end of the options
      if (c == -1)  break;

      switch(c) {
        case 'i':
          std::cout << "-i argument" << optarg << std::endl;
          filenames.push_back(TString(optarg));
          break;
        case 'e':
          n_events = atoi(optarg);
          break;
        case 'h': // Intentional fall-through
        case '?':
          return 0;  
          break;
        default:
          abort();  
      }
    }


    processTChain( filenames, n_events);
    return 0;
}


void processTChain( std::vector<TString> filenames, long n_events ) {

  // create selector instance
  GeneralSkimmer * selector = new GeneralSkimmer();

  // create TChain and add all filenames 
  TChain * tchain = new TChain("Tree");
  for ( std::size_t i = 0; i < filenames.size(); i++ ) {
    tchain->Add( filenames[i] );
  } 

  // process TChain
  tchain->Process(selector, "", n_events);

} 




