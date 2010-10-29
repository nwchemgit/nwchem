// $Id$

/*
 * TokenList is a class for breaking a string into individual tokens
 *           and then analyzing and manipulating them.
 *
 * BGJ (7/00)
 */

#include "TokenList.h"
#include <algorithm>

using namespace std;

// Public member functions

TokenList::TokenList() { }

TokenList::TokenList(const string& Str, const char* Sep)
{
  assign(Str,Sep);
}

/* assign sets the token list to correspond to a given string.  This is
          also useful to be called from the constructor which works with a
          user-supplied string, or from the constructor/assign which works
          with a file after the string has been read. */

void TokenList::assign(const string& Str, const char* Sep)
{
  erase(); // Clear prior contents
  if (Str.length() != 0) {
    OrigString = Str;

    // Make a copy of the original string to break into tokens

    char* TokenString = new char[OrigString.length()+1];
    strcpy(TokenString,OrigString.c_str());

    // Clean up any wackiness before proceeding further

    fixString(TokenString);

    // Get tokens until there are no more

    char* CurToken = strtok(TokenString,Sep);
    while (CurToken != NULL) {
      Tokens.push_back(CurToken);
      CurToken = strtok(NULL,Sep);
    }

    delete [] TokenString;
  }
}

// Assign the token list by reading a line from an input stream

bool TokenList::assign(istream& is, const char* Sep)
{
  bool EoF = 0;
  char Str[512];  // Hopefully this will be long enough!
  if (is.getline(Str,512))
    assign(Str,Sep);
  else {  // Handle EOF case gracefully
    erase();
    EoF = 1;
  }
  return EoF;
}

// A more succint expression of the previous assignment

bool operator>>(istream& is, TokenList& TL)
{
  return TL.assign(is);
}

// Concatentate two TokenLists corresponding to different lines (useful for compound parsing)

void TokenList::operator+=(const TokenList& TL)
{
  copy(TL.Tokens.begin(), TL.Tokens.end(), back_inserter(Tokens));
  // !!! Note that the separator is always a space in this implementation
  OrigString += " " + TL.OrigString;
}

// Query functions

int TokenList::size() const { return Tokens.size(); }
bool TokenList::empty() const { return Tokens.empty(); }
string TokenList::operator[](int i) const { return Tokens[i]; }
string TokenList::getString() const { return OrigString; }

TokenList::const_iterator TokenList::begin() const { return Tokens.begin(); }
TokenList::const_iterator TokenList::end() const { return Tokens.end(); }
TokenList::const_reverse_iterator TokenList::rbegin() const { return Tokens.rbegin(); }
TokenList::const_reverse_iterator TokenList::rend() const { return Tokens.rend(); }

// Private member functions

// Clear out the TokenList

void TokenList::erase()
{
  Tokens.clear();
  OrigString.erase();
}

// Remove any trailing blanks or wacky characters from a string, to give
// the tokenizer a fighting chance

void TokenList::fixString(char* Str)
{
  // Look for the first non-space "normal" character from the end and
  // place '\0' right after it; note that this function does nothing with
  // a null string, which is the correct behavior

  int i;
  for (i = strlen(Str)-1; i >= 0; --i)
    if (isgraph(*(Str+i))) {
      *(Str+i+1) = '\0';
      break;
    }
    else if (i == 0)  // Handle case of all junk
      *Str = '\0';

  // Replace all characters in the "space" class with a literal space,
  // except for newline, which we pass through.  This is done on the
  // assumption that all white space is equivalent.

  for (i = 0; i < strlen(Str); ++i)
    if (isspace(*(Str+i)) && *(Str+i) != '\n')
      *(Str+i) = ' ';
}
