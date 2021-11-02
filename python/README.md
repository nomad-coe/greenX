# pygreenx
*pygreenx* is intended to be a collection of scripts to facilitate the testing and use of the GreenX library.

## Installation
*pygreenx* can be installed from the`greenX/python` directory with:

```bash
pip install -e .
```

and is required to run GreenX's test suite locally. 

# Information For Developers

## External Package Dependencies
If a new external dependency is introduced to the package, this also requires adding to `setup.py` such that pip is aware 
of the new dependency.

## Basic File Structure 
In general, modules should begin with a docstring giving an overview of the module's purpose. External python
libraries should then be imported, followed by a space, then local modules belonging to *pygreenx*.

```angular2html
"""
Functions that operate on lattice vectors 
"""
import numpy as np

from .maths.math_utils import triple_product
```
This may change in the future in favour of loading modules in `__init__.py`. 

## Code Formatting 

We are currently favouring [yapf](https://github.com/google/yapf) formatter, which by default applies PEP8 formatting to 
the code, however any formatter that applies the PEP8 standard is sufficient. 

## Documentation 

### Writing Documentation
All functions and classes should be documented. The favoured docstring is *reStructuredText*:

```angular2html
class SimpleEquation:
   def demo(self, a: int, b: int, c: int) -> list:
    """
    Function definition

    :param int a: quadratic coefficient
    :param int b: linear coefficient 
    :param c: free term
    :type c: int
    :return list y: Function values   
    """
```
where the type can be specified in the `param` description, or separately using the `type` tag. For more details on the
documentation syntax, please refer to this [link](https://devguide.python.org/documenting/).