#ifndef VECOP_H_
#define VECOP_H_

#include "vecop.hpp"

// distance between two vectors
float Distance(vector<float> a, vector<float> b)
{
  if(a.size() == b.size())
  {
    int length = a.size();
    float sum=0;
    for(unsigned int i=0; i<length; i++)
    {   
      sum += pow(a[i] - b[i], 2); 
    }   
    return sqrt(sum);
  }
  else
  {
    VectorError();
  }
}

// angle between three vectors
float Angle(vector<float> a, vector<float> b, vector<float> c)
{
  vector<float> ab, cb;
  ab = VectorSum(a, Scale(b, -1));
  cb = VectorSum(c, Scale(b, -1));

  float ab_cb, ab_ab, cb_cb;
  ab_cb = DotProduct(ab, cb);
  ab_ab = sqrt(DotProduct(ab, ab));
  cb_cb = sqrt(DotProduct(cb, cb));

  float cos_theta = ab_cb / (ab_ab * cb_cb);

  float theta = acos(cos_theta) * 180 / M_PI;

  return theta;
}

// scaling a vector according to a scalar, float
vector<float> Scale(vector<float> a, float mult)
{
  vector<float> a_new;
  for(unsigned int i=0; i<a.size(); i++)
  {
    a_new.push_back(a[i] * mult);
  }
  return a_new;
}

// scale a vector according to a scalar, int
vector<float> Scale(vector<float> a, int mult)
{
  vector<float> a_new;
  for(unsigned int i=0; i<a.size(); i++)
  {
    a_new.push_back(a[i] * mult);
  }
  return a_new;
}

// taking the resultant vector of two vectors
vector<float> VectorSum(vector<float> a, vector<float> b)
{
  vector<float> sum;
  if(a.size() == b.size())
  {
    for(unsigned int i=0; i<a.size(); i++)
    {
      sum.push_back(a[i] + b[i]);
    }
    return sum;
  }
  else
  {
    VectorError();
  }
}

// taking the dot product of two vectors
float DotProduct(vector<float> a, vector<float> b)
{
  float dot_product = 0;
  if(a.size() == b.size())
  {
    for(unsigned int i=0; i<a.size(); i++)
    {
      dot_product += a[i] * b[i];
    }
    return dot_product;
  }
  else
  {
    VectorError();
  }
}

// finding the largest negative integer in an int vector
int FindLargestNegativeInteger(vector<int> a)
{
  int out = a.size() - 1;
  for(unsigned int i=1; i<a.size(); i++)
  {
    if(a[i] > 0)
    {
      out = i - 1;
      break;
    }
  }
  return out;
}

// translating a float vector according to an int scalar
vector<float> Translate(vector<float> start, int offset)
{
  vector<float> end;
  for(unsigned int i=0; i<start.size(); i++)
  {
    end.push_back(start[i] + offset);
  }
  return end;
}

// translating a float vector according to an float scalar
vector<float> Translate(vector<float> start, float offset)
{
  vector<float> end;
  for(unsigned int i=0; i<start.size(); i++)
  {
    end.push_back(start[i] + offset);
  }
  return end;
}

// translating a int vector according to an int scalar
vector<int> Translate(vector<int> start, int offset)
{
  vector<int> end;
  for(unsigned int i=0; i<start.size(); i++)
  {
    end.push_back(start[i] + offset);
  }
  return end;
}

// translating a int vector according to an float scalar
vector<float> Translate(vector<int> start, float offset)
{
  vector<float> end;
  for(unsigned int i=0; i<start.size(); i++)
  {
    end.push_back(start[i] + offset);
  }
  return end;
}

// taking the center of mass of an array of vectors
vector<float> CenterOfMass(vector<vector<float>> coordinates, vector<float> masses)
{
  if(coordinates.size() == masses.size())
  {
    vector<float> center_of_mass;
    float total_mass = VectorSummation(masses);
    for(unsigned int i=0; i<coordinates[0].size(); i++)
    {   
      float sum = 0;
      for(unsigned int j=0; j<coordinates.size(); j++)
      {   
        sum += masses[i] * coordinates[j][i];
      }   
      center_of_mass.push_back(sum/total_mass);
    }   
    return center_of_mass;
  }
  else
  {
    VectorError();
  }
}

// determining whether a value is in an array
bool IsIn(int element, vector<int> array)
{
  for(unsigned int i=0; i<array.size(); i++)
  {
    if(element == array[i])
    {
      return true;
      break;
    }
  }
  return false;
}


// finding a string in a string vector
bool IsIn(string element, vector<string> array)
{
  for(unsigned int i=0; i<array.size(); i++)
  {
    if(element == array[i])
    {
      return true;
      break;
    }
  }
  return false;
}


// find where a string is in a string vector
int Where(string element, vector<string> array)
{
  for(unsigned int i=0; i<array.size(); i++)
  {
    if(element == array[i])
    {
      return i;
      break;
    }
  }
}

// take the sum of all elements in a vector
float VectorSummation(vector<float> a)
{
  float sum=0;
  for(unsigned int i=0; i<a.size(); i++)
  {
    sum += a[i];
  }
  return sum;
}

// vector error message
void VectorError(void)
{
  cout << "Vector Error" << endl;
  exit(1);
}

#endif
