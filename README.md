# CppLightweightMathematics
Simple C++ Vector/Matrix library (useful mainly for graphics-related things)   
Using some state of the art C++ features  

**Status**: unfinished / in development  
**Namespace(s)**: clm  
**Language**: (modern) C++  
  
**Currently offers following functionality**:  
1.  Templated Vector class fully supporting up to 4-dimensional vectors   
2.  Template Matrix class fully supporting up to 4-dimensional matrices  
3.  Additional functions for rad2deg conversion etc.

**Beware of**:  
1.  Matrix & Vector is working with float data and comparing them directly (without epsilon)  
2.  Although some tests were performed, still may contain bugs  

**Usage**:  
Simply include desired headers in your files and you are done.   
All data are in "clm" namespace.   
Both Vector and Matrix class have pre-defined GLSL-like usings by default: vec2/3/4, mat2/3/4 instead of Vector<2/3/4> and Matrix<2/3/4>







