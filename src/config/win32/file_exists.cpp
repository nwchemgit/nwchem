//	$Id$

#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char* argv[])
{
  ifstream is(argv[1]);
  return bool(is);
}
