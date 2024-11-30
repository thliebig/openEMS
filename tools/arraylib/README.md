ArrayLib - Simple Multi-Dimensional Array Library with C-like Syntax
=====================================================================

Motivation
-------------

Originally, multi-dimensional arrays in the simulation engine were
allocated as multi-level C pointers. For example, the following code
allocates and accesses a 3D array:

    // allocate 3D arrays as pointers of pointers
    T*** array=NULL;
    array = new T**[SIZE_I]
    for (size_t i = 0; i < SIZE_I; i++)
    {
    	array[i] = new T*[SIZE_J]
    	for (size_t j = 0; j < SIZE_J; j++)
    	{
            array[i][j] = new T[SIZE_K];
            for (size_t k = 0; k < SIZE_K; k++)
    		{
                // access array via C-style syntax
    			array[i][j][k] = 0;
    		}
    	}
    }

The advantage is that all multi-dimensional arrays can be accessed
like a C-style array. But with several downsides:

1. A single array requires multiple `malloc()` or `new` allocation.
2. The allocated memory is not contiguous, reducing memory locality
and performance, especially when the final dimension is small.

According to the best practice today (e.g. `std::mdspan`), these raw
pointers should be converted to wrapped arrays classes with their
internal indexing code, and they should be accessed using the following
syntax:

    array(i, j, k) = 0;
    array[i, j, k] = 0; // C++23

But this requires massive change of all existing array indexing code
throughout the codebase, which is error-prone and difficult to review.

ArrayLib is a simple library written with the goal of providing a simple
implementation of multi-dimensional array, but preserved the C-like access
syntax using C++ template metaprogramming to minimize the code change in
the simulation engine.

Supported Types
--------------------

1. ArrayIJ - 2D array.
2. ArrayIJK - 3D array.
3. ArrayNIJK - 3D array (N = 3 only), IJK is the coordinate of a 3D
vector in space, N (0, 1, 2) selects the X/Y/Z value from the vector.

To avoid ambiguity, IJK is always used to describe coordinates in
space, not X/Y/Z.

Usage 1: new/delete
--------------------

ArrayLib can be used in several ways, the first method closely mirrors
the C-style code using raw pointers, allowing direct migration from
legacy code.

    ArrayLib::ArrayIJK<float>* arr_ptr = new ArrayLib::ArrayIJK("name", {10, 20, 30});

    // arr_ptr[0][1][2] is wrong, because we need to call the operator[]
    // of ArrayIJK, rather than accessing an array of pointers of type
    // `ArrayIJK*`.
    (*arr_ptr)[0][1][2] = 3;

    // But one can use C++ reference rather than pointers to access it
    // in a readable way, add a line of boilerplate, but more readable
    ArrayLib::ArrayIJK<float>& arr = *arr_ptr;
    arr[0][1][2] = 3;

    delete arr_ptr;

The first method is advantageous in the sense that all existing pointers
with inconsistent lifetimes can be converted to ArrayLib directly as
drop-in replacement. The disadvantage is that the C-style syntax is only
emulated we're accessing a C++ reference, not a C pointer. Yet C++ reference
can't be NULL, so when arrays with inconsistent lifetimes are used in
classes, C-style pointers (with confusing syntax) must be kept as member
variables, not references (with natural syntax).

To write readable code, whenever an array is accessed by a member function,
we must dereference a pointer to create a C++ reference first. Overall, this
is the smoothest migration path, so this style is widely used in engine classes.

Usage 2: As Local Member Variable w/ RAII
------------------------------------------

ArrayLib arrays can be used as local variables by declaring them directly
(with name and size). These arrays are automatically freed by C++'s
RAII mechanics when they go out of scope.

    void process()
    {
        ArrayLib::ArrayIJK<float> arr("name", {10, 20, 30});
        array[0][1][2] = 0;
        // deleted automatically by RAII
    }

This allows one to access arrays using the natural syntax, and without manual
memory management, 

Usage 3: As Class Member Variable w/ RAII
-----------------------------------------

Similarly, ArrayLib arrays can also be used as class member variables
in conformance to the C++ RAII idiom. Class member variables are
recommended for arrays with runtime-determined size.

    class Engine
    {
        // automatically deleteled by RAII
        ArrayLib::ArrayIJK<float> arr1("static size", {10, 20, 30});
        ArrayLib::ArrayIJK<float> arr2;
        
        Engine(size_t i, size_t j, size_t k) :
            arr2("initialize in constructor", {i, j, k})
        {
            // use natural C-style syntax for acess
            arr1[1][2][3] = 0;
            arr2[4][5][6] = 1;
        }
    }

Ideally, this should be the recommended usage for ArrayLib arrays in
C++ classes.

Unfortunately, because the existing codebase does not use the C++ RAII
idiom, and arrays mostly have inconsistent lifetimes (e.g. sometimes,
array sizes can't be determined in the constructor, but instead is
determined in the middle of the class's lifetime from an external
source), it's not possible to translate all existing array usages to
this idiom.

Usage 4: As Class Member Variable w/ Two-Phase Initialization
--------------------------------------------------------------

As a compromise to the lifetime problem, two-phase Initialization is
also supported. This can be done by creating an array object (not
reference or pointers) as a class member variable without specifying
its name and size. The result is a placeholder array in an invalid
state, which can be later initialized.

    class Engine
    {
        // automatically deleteled by RAII, size unknown
        ArrayLib::ArrayIJK<float> arr1;
        
        Engine(size_t i, size_t j, size_t k);
        {
            // two-phase initialization before use
            arr1.Init("name", {i, j, k});

            // use natural C-style syntax for acess
            arr1[1][2][3] = 0;
        }
    }

Although it's potentially unsafe, but one can still enjoy automatic
release of memory upon class destruction by C++'s RAII.

ArrayLib Arrays Can't be Passed by Value
--------------------------------------------

ArrayLib objects are simple wrappers to heap-allocated memory, with few extra
features in comparison to raw C pointers. In fact, they're designed to serve as
drop-in replacement of raw pointers, which is prevalent in the existing codebase
(we expect to refactor the codebase progressively). Thus, they're not meant to
be high-level C++ containers. ArrayLib objects must be always passed as pointers
or references, pass-by-value and copy-assignment are not supported. These
operators are marked as deleted, attempting to do so would create a compile-time
error. They make little sense for the low-level numerical simulation kernel in
its current form.

    ArrayLib::ArrayIJK<float> data("name", size);

    // pass-by-reference, recommended when possible
    void process_ref(ArrayLib::ArrayIJK<float>& arr);
    process_ref(data);

    // pass-by-pointer, use only when necessary
    void process_ptr(ArrayLib::ArrayIJK<float>* arr_ptr);
    process_ptr(&data);

    // pass-by-value, WRONG!
    void process_value(ArrayLib::ArrayIJK<float> arr);
    process_value(data)

Furthermore, practically speaking, in a simulation, large arrays to represent
electromagnetic fields and material properties are allocated and passed around
as pointers. Attempting to pass an array by value is always a bug.
