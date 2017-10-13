/*
 * Copyright: (C) 2014 Walkman Consortium
 * Authors: Juan Alecandro Castano
 * CopyPolicy: Released under the terms of the GNU GPL v2.0.
*/

#ifndef _MATRIX_H2// header guards
#define _MATRIX_H2

//#define M_PI    3.14159
#define MAX_FILTER_LENGTH 8
//typedef float ftype; //number type
typedef double ftype; //number type
#include <cmath>

class FilterH {
	
public:
	//Define Variables
	//************************************************************************************
	int Ncoeff; 					// Number of coefficients
	ftype x[MAX_FILTER_LENGTH]; 	// Circular buffer of incoming values to be filtered
	ftype y[MAX_FILTER_LENGTH];		// Result and accumulators
	ftype a[MAX_FILTER_LENGTH];		// Y coefficient
	ftype b[MAX_FILTER_LENGTH];		// X coefficient
	FilterH();
	/*End Variables definition
	************************************************************************************
	************************************************************************************
	//Begin function heders
	************************************************************************************/
	void clear_filter();  //set values to 0
	
	/*Filter Functions
	T is the sampleTime in second	
	cutoff is the cutoff frequency in Hz
	N is the order of the filter*/

		void least_squares_filter(ftype T, int N); 
		void moving_average_filter(int N);
		void butterworth(ftype T, ftype cutoff, int N);
		void butterDifferentiator(ftype T, ftype cutoff, int N);
        void initialValue(double initialState);

	//Use Filter Function
	ftype applyFilter(double X);
	double differentiator(double X);
	//end Functions
};
#endif
