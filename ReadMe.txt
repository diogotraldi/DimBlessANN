========================================================================
    DYNAMIC LINK LIBRARY : DimBlessANN Project Overview
========================================================================

DimBlessANN.dll is a library that implements an Approximate Nearest Neighbour (ANN)
search algorithm optimized for efficient performance in high-dimensional spaces. 
Unlike traditional methods that struggle in such contexts, this algorithm leverages
the beneficial properties of high dimensions to deliver fast and accurate results.  

This file contains a summary of what you will find in each of the files that
make up your DimBlessANN application.

DimBlessANN.vcxproj
    This is the main project file for VC++ projects generated using an Application Wizard. 
    It contains information about the version of Visual C++ that generated the file, and 
    information about the platforms, configurations, and project features selected with the
    Application Wizard.

DimBlessANN.vcxproj.filters
    This is the filters file for VC++ projects generated using an Application Wizard. 
    It contains information about the association between the files in your project 
    and the filters. This association is used in the IDE to show grouping of files with
    similar extensions under a specific node (for e.g. ".cpp" files are associated with the
    "Source Files" filter).

DimBlessANN.cpp
    This is the main DLL source file.

DimBlessANN.h
    This file contains a class declaration.

AssemblyInfo.cpp
	Contains custom attributes for modifying assembly metadata.

/////////////////////////////////////////////////////////////////////////////
Other notes:

Nearest Neighbour Search is a technique used to find the closest point to a given point in
a multidimensional space. This technique is crucial in various applications such as product
recommendation, pattern recognition, and information retrieval.

DimBlessANN.dll implements an approximate method to perform this search efficiently, 
particularly in high-dimensional spaces where exact methods can be computationally prohibitive.
The algorithm is designed to exploit the favorable characteristics of high dimensions, offering
superior performance in complex scenarios.

/////////////////////////////////////////////////////////////////////////////
