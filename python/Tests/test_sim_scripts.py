# -*- coding: utf-8 -*-
#
# Copyright (C) 2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# The two properties run_testsuite.py relies on in the full-simulation scripts
# next to it. Neither can be seen from a run that passes: a debugging block left
# switched on only blocks the suite on a machine that has a display, and a script
# without a check passes by merely not crashing.

import io
import os
import re
import unittest

TS_DIR = os.path.dirname(os.path.abspath(__file__))


def read(name):
    """The source of one script. Explicitly UTF-8: most of these files carry a
    unit symbol or a Greek letter, and the locale default would be cp1252 or
    worse on a Windows console."""
    with io.open(os.path.join(TS_DIR, name), encoding='utf-8') as fh:
        return fh.read()


def sim_scripts():
    """The full-simulation test scripts of this folder, by name."""
    return sorted(e for e in os.listdir(TS_DIR)
                  if e.endswith('.py')
                  and not e.startswith(('test_', '_'))
                  and e != 'run_testsuite.py')


class Test_SimulationScripts(unittest.TestCase):
    def test_scripts_are_found(self):
        # a typo in the discovery above would make both tests below vacuous
        self.assertTrue(len(sim_scripts()) > 10)

    def test_no_debugging_block_is_switched_on(self):
        # "if 1:" guards the AppCSXCAD viewer and the diagnostic plots. Committed
        # as "if 1:" the viewer waits for a window to be closed, which stalls the
        # whole suite -- and only on a machine with a display, so it survives CI.
        for name in sim_scripts():
            for n, line in enumerate(read(name).splitlines(), 1):
                self.assertIsNone(
                    re.match(r'\s*if 1:', line),
                    '{}:{}: debugging block switched on: {}'.format(
                        name, n, line.strip()))

    def test_every_script_checks_something(self):
        # a script without an assertion cannot fail, so the runner would report
        # PASS for it as long as the simulation does not crash
        for name in sim_scripts():
            self.assertTrue(re.search(r'^\s*assert\s', read(name), re.M),
                            '{}: no assertion, so it can never fail'.format(name))


if __name__ == '__main__':
    unittest.main()
