# MulensModel
Microlensing Modelling package

We present a python package for modelling of microlensing events. The package aims to support a number of higher order effects: microlensing parallax, finite source, binary lens (with lens rotation), and triple lenses in future.

### Second order effects

TBD


### Model conventions

TBD 


### How do we develop the code?

We try to follow [PEP 8](https://www.python.org/dev/peps/pep-0008/) conventions.

How to work on new code?
1. Write use cases if they don't exist yet. 
1. Find the objects, some of them may not be visible in use cases.
1. Organize the objects and how they interact, if needed update documents/Classes.md file.
1. Write a few unit tests _before_ you write any function.
1. Write the function and work on it until unit tests are passed. 


#### Naming conventions
* variables: underscore_separated
* constants: UPPERCASE_AND_UNDERSCORE
* functions: underscore_separated
  * function parameters: underscore_separated
* class names: CamelCaseFirstLetterCapitalized
  * properties: underscore_separated
  * private properties: _leadingUnderscoreAndCamelCase
  * methods: underscore_separated
