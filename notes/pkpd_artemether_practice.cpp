// Step 1: Include required headers and libraries

#include "assert.h" 

/* That line includes the assert.h header, which provides the assert macro in C++.
What assert.h Does: 
    The assert macro is used for debugging by checking assumptions in your code.
    If the condition inside assert(condition) evaluates to false, the program prints an error message and terminates. 
    
When to Use assert:
    Use it to catch programming errors early in development.
    Typically used for checking preconditions in functions.
    It is removed in release builds if NDEBUG is defined (#define NDEBUG).*/


#include "pkpd_artemether.h" // header file for pkpd_artemether class


/* ************************************************************************* */

bool pkpd_artemether::stochastic = true; // setting stochastic to true

/* ************************************************************************* */

// time to declare the constructor
// which means, we are now going to specify what happens when a pkpd_artemether object is created

pkpd_artemether::pkpd_artemether() {
    
}