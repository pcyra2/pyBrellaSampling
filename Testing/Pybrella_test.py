import pyBrellaSampling.pyBrella as pybrella
Test_Dir = "./Testing/TestFiles/pyBrella/"
from unittest.mock import patch
import sys

def test_main():
    args = ["pybrella", "-jt=Test", ]
    with patch.object(sys, "argv", args):
        pybrella.main()

