#include <stdexcept>
#include <string>
using namespace std;

#ifndef RESPONSE_EXCEPTION_H
#define RESPONSE_EXCEPTION_H

class Response_Exception : public logic_error
{
public: 
    Response_Exception(const string& message = "") : logic_error(message.c_str()) {}	
};

#endif
