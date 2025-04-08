# pkpd_artemether:

## Required files and corresponding descriptions:

### main.cpp: 
	Contains the main() function, where the program runs and interacts with the class.

### pkpd_artemether.cpp: 
	Contains the implementation of the functions declared in the class (pkpd_artemether).

### pkpd_artemether.h: 
	Declares the class interface (functions and variables) without implementation details.

<br>

## pkpd_artemether.cpp structure/control of flow:

1. Load required headers and libraries
	-
	* ``` #include "assert.h" ``` : 

		* Provides the assert macro in C++.
		* The assert macro is used for debugging by checking assumptions in your code.
		* If the condition inside assert(condition) evaluates to false, the program prints an error message and terminates. 

    	* Use it to catch programming errors early in development.
		* Typically used for checking preconditions in functions.
    	* It is removed in release builds if NDEBUG is defined (#define NDEBUG).

	<br>

	* ``` #include "pkpd_artemether.h" ```: 

		* Header file outlining the pkpd_artemether class along with required variables and functions.

		* We only declare the variables in the header file.

<br>

2. Set bool ```pkpd_artemether::stochastic``` to ```true```
	- 
<br>

3. Define the constructor function for ```pkpd_artemether```
	-
	* Name of the constructor function: ```pkpd_artemether()```

		* Constructor function **declared** in ```pkpd_artemether.h```

		* Constructor **declaration** tells the compiler that a constructor exists, but does not specify what it does.
	
	* Constructor **definition** implements what the constructor does.
		* Definition specifies what happens when an object of pkpd_artemether is created.
		
