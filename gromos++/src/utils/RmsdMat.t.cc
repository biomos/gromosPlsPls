#include <iostream>
#include "RmsdMat.h"
#define MATSIZ 100


int main(){

  unsigned int i, j;
  utils::RmsdMat mat1(MATSIZ), mat2(MATSIZ);

  
  cout << "inserting" << endl;
  for(i = 0; i < MATSIZ; i++){
    for(j = i + 1; j < MATSIZ; j++){
      cout << i << "\t" << j << "\t" << i + j << endl;
      mat1.insert(i, j, float(i + j));
    }
  }

  cout << "transposing" << endl;
  // make mat2 the transposed of mat1
  for(i = 0; i < MATSIZ; i++){
    for(j = i + 1; j < MATSIZ; j++){
      cout << i << "\t" << j << "\t" << i + j << endl;
      mat2.insert(i, j, mat1.retrieve(i, j));
    }
  }

  cout << "extracting" << endl;
  // compare (should be identical)
  for(i = 0; i < MATSIZ; i++){
    for(j = i + 1; j < MATSIZ; j++){
      if(mat1.retrieve(i, j) != mat2.retrieve(i, j)){
        exit(1);
      }
    }
  }
  return 0;
}
