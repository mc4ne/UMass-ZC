#ifndef UMASS_ALGORITHM_HH
#define UMASS_ALGORITHM_HH

namespace UMass {
   namespace Algorithm { 
      //______________________________________________________________________________
      template <typename T>
	 void BinarySearch(std::vector<T> array,T key,int &lowerbound,int &upperbound){
	    const int N         = array.size();
	    int comparisonCount = 1;    //count the number of comparisons (optional)
	    lowerbound          = 0;
	    upperbound          = N-1;  
	    // To start, find the subscript of the middle position.
	    int position = (lowerbound + upperbound) / 2;
	    while( (array[position]!=key) && (lowerbound <= upperbound)){
	       comparisonCount++;
	       if (array[position] > key){
		  // decrease position by one.
		  upperbound = position - 1;
	       }else{
		  // else, increase position by one.
		  lowerbound = position + 1;
	       }  
	       position = (lowerbound + upperbound) / 2;
	    }  
	    // T lo=0,hi=0,mid=0;
	    int dump = lowerbound;
	    if (lowerbound <= upperbound){
	       // Here we have an exact match to the key
	       lowerbound = position;
	       upperbound = position+1;
	    }else{
	       // Here the lower bound surpassed the upper
	       lowerbound = upperbound;
	       upperbound = dump;
	    }  
	    // to safeguard against values that are outside the boundaries of the grid 
	    if(upperbound>=N){
	       upperbound = N-1;
	       lowerbound = N-2;
	    }  
	    if(upperbound==0){
	       lowerbound = 0;
	       upperbound = 1;
	    }
	 }
   } //::Algorithm 
} //::UMass

#endif 
