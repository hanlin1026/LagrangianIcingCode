#include <findAll.h>

using namespace std;

void findAll(vector<int>& One, vector<int>& Match, int Num) {
  // Function to compare two vectors of indices and return a 
  // vector listing the indices of One where common elements exist
  
  for (int i=0; i<One.size(); i++) {
    if (One[i] == Num) {
      Match.push_back(i);
    }
  }
}
