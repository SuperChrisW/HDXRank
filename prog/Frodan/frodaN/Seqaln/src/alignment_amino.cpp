///A tutorial about heuristic multiple sequence alignment.
#include <iostream>
#include <seqan/align.h>
#include <seqan/graph_msa.h>
#include <seqan/file.h>
#include <string>
#include <cstring>
#include <cstdlib>
using namespace seqan;
using namespace std;

int main(int argc, const char *argv[])
{
  string line;
  string seq[2];
  int i=-1;
  ifstream fstrm(argv[1], ios_base::in | ios_base::binary);
  
  while(!fstrm.eof()){
    getline( fstrm, line);
    if ( line[0] != '>'){
      seq[i] +=   line;

    }
    else{
      i++;
    }
  }  
  fstrm.close();
  if ( i<1 )
    cout << "File "<<argv[1]<<" does not have two sequences in Fasta format"<<endl;

/// A set of sequences to align
  typedef String<AminoAcid> TSequence;

  TSequence seq1= seq[0];
  TSequence seq2= seq[1];

  /*
   //If sequences are in 2 separate files uncomment this lines
  ::std::fstream fstrm1;
  fstrm1.open(argv[1], ::std::ios_base::in | ::std::ios_base::binary);
  read(fstrm1, seq1, Fasta());
  
  fstrm1.close();
  ::std::fstream fstrm2;
  fstrm2.open(argv[2], ::std::ios_base::in | ::std::ios_base::binary);
  read(fstrm2, seq2, Fasta());
  fstrm2.close();
  */
  
/// Scoring: Blosum62 where gex = -1, gop = -11
    Blosum62 sc(-1, -11);
    
    /// Create an Align object or alternatively an Alignment graph object to store the MSA
      Align<TSequence, ArrayGaps> align;
      resize(rows(align), 2);
      assignSource(row(align, 0), seq1);
      assignSource(row(align, 1), seq2);
      //::std::cout << align<< ::std::endl;
      //
      
      /// Heuristic PSA
	int score = globalAlignment(align, sc);
		
      /// Output of the PSA
	::std::cout <<"Score: " <<score << ::std::endl;
      //::std::cout << align << ::std::endl;
	
	::std::ofstream fstrm_out(argv[2],::std::ios::out | ::std::ios::binary);
	fstrm_out << align;
	
	return 0;
}
