#pragma once
#include <string>

class LibpllException: public std::exception {
public:
  explicit LibpllException(const std::string &s): msg_(s) {}
  explicit LibpllException(const std::string &s1, 
      const std::string &s2): msg_(s1 + s2) {}
  virtual const char* what() const noexcept { return msg_.c_str(); }
  virtual ~LibpllException() {}
private:
  std::string msg_;
};

