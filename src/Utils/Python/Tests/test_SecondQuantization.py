__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

import pytest
import scine_utilities as scine
import scine_utilities.sq as sq
import numpy as np
import os

def test_ElectronicDeterminant():
    det = sq.ElectronicDeterminant("2a00")
    assert det.is_occupied(sq.SpinComponent.Alpha, 0) == True
    assert det.is_occupied(sq.SpinComponent.Beta, 0) == True
    assert det.is_occupied(sq.SpinComponent.Alpha, 1) == True
    assert det.is_occupied(sq.SpinComponent.Beta, 1) == False
    assert det.is_occupied(sq.SpinComponent.Alpha, 2) == False
    assert det.is_occupied(sq.SpinComponent.Beta, 2) == False
    assert det.is_occupied(sq.SpinComponent.Alpha, 3) == False
    assert det.is_occupied(sq.SpinComponent.Beta, 3) == False

    result = det.count_electrons()
    assert result == (2,1)

    assert det.occupied(sq.SpinComponent.Alpha) == [0,1]
    assert det.occupied(sq.SpinComponent.Beta) == [0]
    assert det.virtual(sq.SpinComponent.Alpha) == [2,3]
    assert det.virtual(sq.SpinComponent.Beta) == [1,2,3]

def test_ExciteElectronicDeterminant():
    det = sq.ElectronicDeterminant("2200")
    exc = sq.SingleExcitation(1,3,sq.SpinComponent.Alpha)
    det.apply_excitation(exc.occ, exc.vir, exc.spin)
    assert det.is_occupied(sq.SpinComponent.Alpha, 1) == False
    assert det.is_occupied(sq.SpinComponent.Alpha, 3) == True

    det.apply_cross_excitation(3,3,sq.SpinComponent.Alpha, sq.SpinComponent.Beta)
    assert det.is_occupied(sq.SpinComponent.Alpha, 3) == False
    assert det.is_occupied(sq.SpinComponent.Beta, 3) == True

def test_FlipOnv():
    onv = sq.OccupationNumberVector([1,1,0,0])
    onv.flip(3)
    assert onv.count_occupied() == 3
    assert onv.count_virtual() == 1
    assert onv.size() == 4
    assert onv.is_occupied(3) == True
