"""
Optimnus binding unit testing
-----------------------------

Version: 0.0.1

Authors:

etc. etc. 

Using pytest to do the following

- Test main binaries (foldx, ialign) 

- Test mutation file 

For this purpose, I have copied 5 pdb files as input in this folder as follows: [] 


"""

import os
import sys
import pytest

sys.path.append(os.path.abspath(os.path.join('..', 'src'))) # Reference to the src path

from mutation import foldx_align # python code containing the foldx mutation function

@pytest.fixture()
def test_binaries():
	pass
	

def test_callfoldx():
	pass


	
