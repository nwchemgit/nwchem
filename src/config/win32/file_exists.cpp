//	$Id: file_exists.cpp,v 1.1 2000-07-27 15:54:44 bjohnson Exp $

#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char* argv[])
{
  ifstream is(argv[1]);
  return bool(is);
}
