#include <vector>
#include <iostream>
#include <cmath>
#include <string>

using namespace std;

// distance between two vectors
float Distance(vector<float>, vector<float>);

// angle between three vectors
float Angle(vector<float>, vector<float>, vector<float>);

// scaling a vector according to a scalar
vector<float> Scale(vector<float>, float);
vector<float> Scale(vector<float>, int);

// taking the resultant vector of two vectors
vector<float> VectorSum(vector<float>, vector<float>);

// taking the dot product of two vectors
float DotProduct(vector<float>, vector<float>);

// finding the largest negative integer in an int vector
int FindLargestNegativeInteger(vector<int>);

// translating a vector according to a scalar
vector<float> Translate(vector<float>, int);
vector<float> Translate(vector<float>, float);
vector<int> Translate(vector<int>, int);
vector<float> Translate(vector<int>, float);

// taking the center of mass of an array of vectors
vector<float> CenterOfMass(vector<vector<float>>, vector<float>);

// determining whether a value is in an array
bool IsIn(int, vector<int>);
bool IsIn(string, vector<string>);

// finding the index of a value in an array
int Where(string, vector<string>);

// taking the summation of values in a vector
float VectorSummation(vector<float>);

// printing vector error
void VectorError(void);
