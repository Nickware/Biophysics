"""Some definitions."""

"""Definitions for my examples."""

import os, os.path

from pCore import TestScript_InputDataPath, \
                  TestScript_OutputDataPath

_name = "myExamples"

dataPath = TestScript_InputDataPath(_name)
outPath  = TestScript_OutputDataPath(_name)

xyzPath = os.path.join(dataPath, "xyz")

_FullVerificationSummary = False
