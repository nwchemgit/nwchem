#ifndef _H_TOKENLIST
#define _H_TOKENLIST

// $Id$

/*
 * TokenList is a class for breaking a string into individual tokens
 *           and then analyzing and manipulating them.
 *
 * BGJ (7/00)
 */

#pragma warning( disable : 4786 )

#include <string>
#include <iosfwd>
#include <vector>

class TokenList {
public:
  typedef std::vector<std::string>::const_iterator const_iterator;
  typedef std::vector<std::string>::const_reverse_iterator const_reverse_iterator;

  TokenList();
  TokenList(const std::string& Str, const char* Sep=" ");

  void assign(const std::string& Str, const char* Sep=" ");
  bool assign(std::istream& is, const char* Sep=" ");
  friend bool operator>>(std::istream& is, TokenList& TL);

  // Concatentate two TokenLists corresponding to different lines (useful for compound parsing)

  void operator+=(const TokenList& TL);

  // Query functions

  int size() const;
  bool empty() const;
  std::string operator[](int i) const;
  std::string getString() const;

  // Iterator functions

  const_iterator begin() const;
  const_iterator end() const;
  const_reverse_iterator rbegin() const;
  const_reverse_iterator rend() const;

private:
  std::string OrigString;
  std::vector<std::string> Tokens;

  void erase();
  static void fixString(char* Str);
};

#endif /* _H_TOKENLIST */
