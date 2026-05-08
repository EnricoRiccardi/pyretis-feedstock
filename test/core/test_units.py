# Copyright (c) 2026, PyRETIS Development Team.
# Distributed under the LGPLv2.1+ License. See LICENSE for more info.
"""Test the functionality of units from pyretis.core.

We also do some test to check the correctness of values.
"""
from copy import deepcopy
import logging
import os
import tempfile
import pytest
from copy import deepcopy
import logging
import os
import tempfile
from pyretis.core.units import (
    create_conversion_factors,
    generate_system_conversions,
    generate_conversion_factors,
    read_conversions,
    CONVERT,
    UNITS,
    write_conversions,
    units_from_settings,
    convert_bases,
    bfs_convert,
)
logging.disable(logging.CRITICAL)


def test_create_lennard_jones_units(caplog):
    """Test that we create correct Lennard-Jones units."""
    create_conversion_factors('lj', length=(3.405, 'A'),
                              energy=(119.8, 'kB'),
                              mass=(39.948, 'g/mol'),
                              charge='e')
    assert CONVERT['length']['bohr', 'nm'] == pytest.approx(
        0.052917721067000, abs=1e-12
    )
    assert CONVERT['length']['lj', 'nm'] == pytest.approx(0.34050, abs=1e-4)
    assert CONVERT['length']['bohr', 'lj'] == pytest.approx(
        0.155411809301, abs=1e-12
    )
    assert CONVERT['time']['lj', 'ps'] == pytest.approx(
        2.156349772323142, abs=1e-12
    )
    assert CONVERT['mass']['lj', 'kg'] == pytest.approx(
        6.633521358698435e-26, abs=1e-12
    )
    assert CONVERT['force']['lj', 'N'] == pytest.approx(
        4.857612120293684e-12, abs=1e-12
    )
    assert CONVERT['energy']['lj', 'J'] == pytest.approx(
        1.654016926960000e-21, abs=1e-12
    )
    assert CONVERT['velocity']['lj', 'nm/ps'] == pytest.approx(
        1.579057369867981e-01, abs=1e-12
    )
    assert CONVERT['pressure']['lj', 'bar'] == pytest.approx(
        4.189754740302599e+02, abs=1e-12
    )
    assert CONVERT['temperature']['lj', 'K'] == pytest.approx(
        1.198000000000000e+02, abs=1e-12
    )
    assert CONVERT['charge']['lj', 'e'] == pytest.approx(4.940801883873017e-02,
                                                         abs=1e-12)


def test_create_cgs_units(caplog):
    """Test that we create correct cgs units."""
    create_conversion_factors('cgs', length=(0.01, 'm'),
                              energy=(1.0e-7, 'J'),
                              mass=(1.0, 'g'), charge='e')
    assert CONVERT['force']['cgs', 'dyn'] == pytest.approx(1.0, abs=1e-12)
    assert CONVERT['force']['cgs', 'N'] == pytest.approx(1.0e-5, abs=1e-12)
    assert CONVERT['time']['cgs', 's'] == pytest.approx(1.0, abs=1e-12)
    assert CONVERT['length']['cgs', 'm'] == pytest.approx(1.0e-2, abs=1e-12)
    assert CONVERT['mass']['cgs', 'kg'] == pytest.approx(1.0e-3, abs=1e-12)
    assert CONVERT['velocity']['cgs', 'm/s'] == pytest.approx(1.0e-2,
                                                              abs=1e-12)
    assert CONVERT['energy']['cgs', 'J'] == pytest.approx(1.0e-7, abs=1e-12)
    assert CONVERT['pressure']['cgs', 'Pa'] == pytest.approx(1.0e-1, abs=1e-12)
    assert CONVERT['temperature']['cgs', 'K'] == pytest.approx(1.0, abs=1e-12)
    assert CONVERT['charge']['cgs', 'C'] == pytest.approx(
        3.335640951981520e-10, abs=1e-12
    )


def test_have_common_unit_systems(caplog):
    """Test that we can generate all conversions."""
    systems = ('lj', 'real', 'metal', 'au',
               'electron', 'si', 'gromacs')
    for sys in systems:
        create_conversion_factors(sys)
    all_pairs = []
    for sys1 in systems:
        for sys2 in systems:
            if sys1 != sys2:
                generate_system_conversions(sys1, sys2)
                all_pairs.append((sys1, sys2))
    for key in CONVERT:
        msg = ['Could not find conversion "{}" -> "{}"',
               f'for dimension "{key}"']
        for pair in all_pairs:
            msg[0] = msg[0].format(*pair)
            msgtxt = ' '.join(msg)
            assert pair in CONVERT[key], msgtxt


def test_creation_of_units(caplog):
    """Test that creation of units works and fails as expected."""
    with pytest.raises(ValueError):
        create_conversion_factors('test', length=None, energy=None,
                                  mass=None, charge=None)
    with pytest.raises(ValueError):
        create_conversion_factors('test', length=(1.0, 'm'),
                                  energy=(1.0, 'J'), mass=(1.0, 'kg'),
                                  charge=None)
    with pytest.raises(ValueError):
        create_conversion_factors('test', length=(1.0, 'm'),
                                  energy=(1.0, 'J'), mass=(1.0, 'kg'),
                                  charge='a non-existing unit')
    with pytest.raises(LookupError):
        create_conversion_factors('test', length=(1.0, 'strange_unit'),
                                  energy=(1.0, 'J'), mass=(1.0, 'kg'),
                                  charge='e')
    with pytest.raises(TypeError):
        create_conversion_factors('test', length=1.0,
                                  energy=(1.0, 'J'), mass=(1.0, 'kg'),
                                  charge='e')
    # the next one should be successful
    create_conversion_factors('test', length=(1.0, 'm'),
                              energy=(1.0, 'J'), mass=(1.0, 'kg'),
                              charge='e')
    generate_system_conversions('test', 'real')
    # check if we indeed created all conversions
    for key in CONVERT:
        dimtxt = f'for dimension "{key}"'
        for unit in UNITS[key]:
            pair = ('test', unit)
            msg = 'Could not find conversion "{}" -> "{}"'.format(*pair)
            msgtxt = ' '.join([msg, dimtxt])
            assert pair in CONVERT[key], msgtxt
            pair = (unit, 'test')
            msg = 'Could not find conversion "{}" -> "{}"'.format(*pair)
            msgtxt = ' '.join([msg, dimtxt])
            assert pair in CONVERT[key], msgtxt

    with pytest.raises(ValueError):
        create_conversion_factors('test', length=(1.0, 'm'),
                                  energy=(1.0, 'J'),
                                  mass=(1.0, 'kg'), charge=(100, 'e'))


def test_read_from_file(caplog):
    """Test that we can read units from a input file."""
    dirname = os.path.dirname(os.path.abspath(__file__))
    filename = os.path.join(dirname, 'units_input.txt')
    conv = read_conversions(filename=filename,
                            select_units='test_system')
    assert conv['charge']['C', 'test_system'] == pytest.approx(
        1.23450, abs=1e-12
    )
    assert conv['energy']['J', 'test_system'] == pytest.approx(
        6.7890, abs=1e-12
    )
    assert conv['force']['N', 'test_system'] == pytest.approx(
        4.2424242e1, abs=1e-12
    )
    assert conv['length']['m', 'test_system'] == pytest.approx(
        1.111112e1, abs=1e-12
    )
    assert conv['mass']['kg', 'test_system'] == pytest.approx(
        1.234e-4, abs=1e-12
    )
    assert conv['pressure']['Pa', 'test_system'] == pytest.approx(
        1e-1, abs=1e-12
    )
    assert conv['temperature']['K', 'test_system'] == pytest.approx(
        1.0, abs=1e-12
    )
    assert conv['time']['s', 'test_system'] == pytest.approx(0.5, abs=1e-12)
    assert conv['velocity']['m/s', 'test_system'] == pytest.approx(
        2.0, abs=1e-12
    )


def test_write_conversions(caplog):
    """Test that we can write out the conversions."""
    with tempfile.NamedTemporaryFile() as temp:
        write_conversions(filename=temp.name)
        temp.flush()
        conv = read_conversions(filename=temp.name)
        for key, val in conv.items():
            assert key in CONVERT
            assert val == pytest.approx(CONVERT[key])


def test_convert_bases(caplog):
    """Test convert_bases function"""
    copy_convert = deepcopy(CONVERT)
    CONVERT['mass'] = {('kg', 'kg'): 1.0, ('g', 'kg'): 0.001,
                       ('g', 'g'): 1.0}
    with caplog.at_level(logging.WARNING):
        convert_bases('mass')
    assert CONVERT['mass'][('kg', 'g')] == pytest.approx(1000.0)
    CONVERT.clear()
    CONVERT.update(copy_convert)


def test_bfs_convert(caplog):
    """Test bfs_convert function"""
    convertdim = {('g', 'g'): 1.0, ('kg', 'g'): 1000.0, ('kg', 'kg'): 1.0}
    assert bfs_convert(convertdim, 'g', 'kg') == (('g', 'kg'), None, None)
    assert bfs_convert(convertdim, 'kg', 'kg') == (('kg', 'kg'), 1.0, None)
    assert bfs_convert(convertdim, 'kg', 'g') == (('kg', 'g'), 1000.0, None)


def test_generate_conversion(caplog):
    """Test the generator of conversion factors."""
    generate_conversion_factors('myunit', (1.0, 'm'), (1.0, 'J'),
                                (1.0, 'kg'), charge='e')
    assert CONVERT['mass']['myunit', 'g'] == pytest.approx(1000.)
    generate_conversion_factors('myunit2', (1.0, 'm'), (1.0, 'J'),
                                (1.0, 'kg'), charge='C')
    # Test with missing kB
    with pytest.raises(ValueError):
        generate_conversion_factors('myunit3', (1.0, 'm'), (1.0, 'JJ'),
                                    (1.0, 'kg'), charge='e')
    # Test with non-unit "1.0" values
    generate_conversion_factors('myunit4', (2.0, 'm'), (0.5, 'J'),
                                (3.0, 'g'))
    assert CONVERT['length']['myunit4', 'm'] == pytest.approx(2.0)
    assert CONVERT['length']['m', 'myunit4'] == pytest.approx(0.5)
    assert CONVERT['energy']['myunit4', 'J'] == pytest.approx(0.5)
    assert CONVERT['energy']['J', 'myunit4'] == pytest.approx(2.0)
    assert CONVERT['mass']['myunit4', 'g'] == pytest.approx(3.0)
    assert CONVERT['mass']['g', 'myunit4'] == pytest.approx(0.33333333)
    assert CONVERT['mass']['myunit4', 'kg'] == pytest.approx(0.003)
    # Test with 'unit' label part of  UNITS[dim]. To check coverage.
    copy_convert = deepcopy(CONVERT)
    assert CONVERT == copy_convert
    generate_conversion_factors('bohr', (1., 'm'), (1., 'J'), (1., 'g'))
    assert CONVERT != copy_convert
    assert CONVERT['length']['bohr', 'bohr'] == pytest.approx(1.)
    # Test some expected but wrong assignment
    assert CONVERT['length']['bohr', 'm'] == pytest.approx(1.)
    assert CONVERT['mass']['bohr', 'g'] == pytest.approx(1.)
    # Restore CONVERT
    CONVERT.clear()
    CONVERT.update(copy_convert)
    assert CONVERT == copy_convert
    with pytest.raises(Exception) as context:
        generate_conversion_factors('xx_au', (1., 'm'), (1., 'J'), (1., 'xx'))
    assert "Missing" in str(context.value)


def test_units_from_settings(caplog):
    """Test that we can create units from settings."""
    # Using a predefined unit system
    settings = {'system': {'units': 'real'}}
    msg = units_from_settings(settings)
    assert msg == 'Created units: "real".'
    # Using a custom system:
    settings = {'system': {'units': 'my_new_system1'}}
    settings['unit-system'] = {
        'name': 'my_new_system1',
        'length': (1.0, 'bohr'),
        'mass': (1.0, 'g'),
        'energy': (1.0, 'J'),
        'charge': 'e'
    }
    msg = units_from_settings(settings)
    assert msg == 'Created unit system: "my_new_system1".'
    # Test for errors:
    with pytest.raises(ValueError):
        settings = {'system': {'units': 'my_new_system2'}}
        settings['unit-system'] = {
            'length': (1.0, 'bohr'),
            'mass': (1.0, 'g'),
            'energy': (1.0, 'J'),
            'charge': 'e'
        }
        units_from_settings(settings)
    # Inconsistent name:
    with pytest.raises(ValueError):
        settings = {'system': {'units': 'my_new_system2'}}
        settings['unit-system'] = {
            'name': 'my_new_system3',
            'length': (1.0, 'bohr'),
            'mass': (1.0, 'g'),
            'energy': (1.0, 'J'),
            'charge': 'e'
        }
        units_from_settings(settings)
    # Missing one of length, mass, energy, charge
    for i in ('length', 'mass', 'energy', 'charge'):
        settings = {'system': {'units': 'my_new_system2'}}
        settings['unit-system'] = {
            'name': 'my_new_system2',
            'length': (1.0, 'bohr'),
            'mass': (1.0, 'g'),
            'energy': (1.0, 'J'),
            'charge': 'e'
        }
        del settings['unit-system'][i]
        with pytest.raises(ValueError):
            units_from_settings(settings)
