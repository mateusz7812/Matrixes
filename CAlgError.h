#pragma once
#include <vcruntime_exception.h>

class BAD_SIZE_EXCEPTION : public std::exception
{
    const char* what() const throw ()
    {
        return "bad size exception";
    }
};

class BAD_VALUE_EXCEPTION : public std::exception
{
    const char* what() const throw ()
    {
        return "bad value exception";
    }
};

class SINGULAR_MATRIX_EXCEPTION : public std::exception
{
    const char* what() const throw ()
    {
        return "singular matrix exception";
    }
};

