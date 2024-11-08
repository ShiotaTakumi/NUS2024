#! /usr/bin/env python

import numpy as np

# Check if float128 is supported
try:
    x = np.float128(1.0)
    print(f"float128 is available: {x}\n")
except AttributeError:
    print("float128 is not available.")
