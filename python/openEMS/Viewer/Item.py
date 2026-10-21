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

"""openEMS.Viewer.Item
Common bases of the field items: :class:`FieldItem` (colour map and legend,
linear/log scale, value range, tree groups and drawing units) and
:class:`DumpItem` (a field dump: quantity, frequency or time step, phase or
time animation and normalisation).
"""

import numpy as np

from vtkmodules.vtkCommonCore import vtkLookupTable

from CSXCAD.Viewer import ViewItem

from openEMS.Viewer import Colors, Quantity

#: labels of the result groups in the scene tree
GROUP_LABELS = {'results': 'Results', 'results/fields': 'Fields', 'results/farfield': 'Far field'}


class FieldItem(ViewItem):
    """ Base of the field items, see :class:`FieldSlices`. """

    HELP = dict(
        scale='Linear values, or dB: 20 log10 of field values, 10 log10 of power values such as SAR.',
        dynamic_range='dB scale: how far below the maximum the colours reach; lower values get the lowest colour.',
        value_range='Colour range: "auto" (from the data; symmetric for signed values, the same over all\n'
                    'animation steps) or two numbers "min, max".',
        cmap='Colour map (a matplotlib name).',
        legend='Show the colour legend. Drag it with the mouse; a double click puts it back.',
    )

    def __init__(self, name, cmap=None, scale='linear', dynamic_range=40.0, value_range='auto', legend=True):
        super().__init__(name, include_in_bounds=True)
        self.legend = bool(legend)
        self._lut = vtkLookupTable()
        self._bar = Colors.ScalarBar(self._lut)
        self.AddProp(self._bar)
        self.cmap = cmap
        self.scale = self._CheckScale(scale)
        self.dynamic_range = float(dynamic_range)
        self.value_range = Quantity.ParseRange(value_range)
        self._delta = None    # drawing unit in m, from the viewer's structure

    @staticmethod
    def _CheckScale(scale):
        if scale not in ('linear', 'log'):
            raise ValueError('scale must be "linear" or "log"')
        return scale

    # -- viewer hooks
    def OnAdded(self, viewer):
        group = self.GetGroup()
        for path, label in GROUP_LABELS.items():
            if group == path or group.startswith(path + '/'):
                viewer.AddGroup(path, label)
        viewer.AddLegend(self, self._bar, render=False)
        csx = viewer.GetCSX()
        if csx is not None:
            self._delta = csx.GetGrid().GetDeltaUnit()
        self.Update()

    def Update(self):
        """ Rebuild the props from the settings (must not render). """

    # -- shared settings
    def _Cmap(self, signed):
        return self.cmap or Colors.DefaultCmap(self._Kind(), signed)

    def _Kind(self):
        return ''

    def _ApplyLut(self, signed, vrange):
        Colors.LookupTable(self._Cmap(signed), lut=self._lut)
        self._lut.SetTableRange(*vrange)
        self._ShowLegend()

    def _ShowLegend(self, wanted=True):
        """ Show the colour legend if the setting, `wanted` and the item say so. """
        self._bar.SetVisibility(self.legend and wanted and self.GetVisible())

    def _CommonParameters(self, signed, scale_default='linear'):
        cmaps = list(Colors.CMAPS)
        cmap = self._Cmap(signed)
        if cmap not in cmaps:
            cmaps.append(cmap)
        vr = self.value_range
        return [
            dict(name='scale', type='choice', value=self.scale, choices=['linear', 'log'], default=scale_default,
                 section='Colours'),
            dict(name='dynamic_range', label='dynamic range (dB)', type='float', value=self.dynamic_range,
                 min=1.0, max=200.0, step=5.0, default=40.0, section='Colours'),
            dict(name='value_range', label='value range', type='str', value=vr, default='auto',
                 section='Colours'),
            dict(name='cmap', label='colour map', type='choice', value=self.cmap or cmap, choices=cmaps,
                 default=self.cmap or cmap, section='Colours'),
            dict(name='legend', label='colour legend', type='bool', value=self.legend, default=True,
                 section='Colours'),
        ]

    def _SetCommon(self, name, value):
        """ Set a shared setting; False if `name` is none. """
        if name == 'scale':
            self.scale = self._CheckScale(value)
        elif name == 'dynamic_range':
            self.dynamic_range = float(value)
        elif name == 'value_range':
            self.value_range = Quantity.ParseRange(value)
        elif name == 'cmap':
            Colors.LookupTable(value)   # validates
            self.cmap = value
        elif name == 'legend':
            self.legend = bool(value)
        else:
            return False
        return True

    def _Range(self, peak, signed):
        if self.value_range != 'auto':
            return tuple(self.value_range)
        return Quantity.AutoRange(peak, signed, self.scale, Quantity.IsPower(self._Kind()), self.dynamic_range)

    def _LengthText(self, pos):
        unit = Quantity.LengthUnit(self._delta or 1.0)
        return '({}){}'.format(', '.join('{:.3g}'.format(c) for c in pos), ' ' + unit if unit else '')


def ParseNorm(value):
    """ A normalisation factor from a number or its text, e.g. ``'2-1j'``. """
    if isinstance(value, str):
        value = complex(value.replace(' ', ''))
    value = complex(value)
    if value == 0:
        raise ValueError('the normalisation factor must not be 0')
    return value.real if value.imag == 0 else value


class DumpItem(FieldItem):
    """ Base of the items showing a field dump (:class:`FieldSlices`,
    :class:`FieldIsoSurface`): the source and quantity, the shown frequency
    or time step, the phase or time animation and the normalisation. A
    subclass reads its data through :meth:`_Norm` and implements
    :meth:`_ShowStep`, :meth:`_TimePeaks` and :meth:`_SetOwn`. """

    #: default of `amplitude`: None follows the quantity (the magnitude and
    #: scalars are shown as amplitude, components as instantaneous value)
    AMPLITUDE_DEFAULT = None
    #: settings that change the read data
    DATA_SETTINGS = ('quantity', 'domain', 'freq', 'sample', 'time_sample', 'norm')

    HELP = dict(
        FieldItem.HELP,
        quantity='What is shown: a component (x, y, z; rho, alpha on a cylindrical mesh), the magnitude,\n'
                 'or the value of a scalar dump such as SAR.',
        domain='Frequency domain (at one frequency; a time-domain dump is transformed) or time domain\n'
               '(the recorded time steps).',
        freq='The frequency shown: a stored one, or any frequency for a time-domain dump.',
        phase='Phase in degrees of the instantaneous value Re(F exp(j phase)); the animation starts there.',
        amplitude='Show the amplitude |F| instead of the instantaneous value at the phase.',
        time_step='The recorded time step shown. It follows the animation, and setting it moves the animation.',
        sample='Use every n-th mesh line only: less memory and faster for large dumps.',
        norm='The values are divided by this factor, e.g. the incident port voltage at the frequency.\n'
             'A complex factor such as 2-1j also removes its phase.',
        unit='Unit of the normalised values, e.g. V/m/V; empty: the unit of the dump.',
        animate='Animate over the phase (frequency domain) or the time steps (time domain).',
        n_phase='Animation steps per period of the phase.',
        time_sample='Animate every n-th recorded time step.',
    )

    def __init__(self, name, source, quantity=None, domain=None, freq=None, time_step=None, phase=0.0,
                 amplitude=None, sample=1, time_sample=1, n_phase=36, animate=True, norm=1.0, unit=None,
                 cmap=None, scale='linear', dynamic_range=40.0, value_range='auto', legend=True):
        super().__init__(name, cmap, scale, dynamic_range, value_range, legend)
        self.source = src = source
        self._choices = Quantity.Components(src.mesh_type, src.is_vector)
        self.quantity = Quantity.Normalize(quantity, src.kind, self._choices)
        self.domain = domain or ('frequency' if src.HasFD() else 'time')
        self._CheckDomain(self.domain)
        self.freq = self._SnapFreq(freq) if freq is not None else (float(src.freqs[0]) if src.HasFD() else None)
        if self.domain == 'frequency' and self.freq is None:
            raise ValueError('give `freq` to show a time-domain dump at a frequency')
        self.time_step = None if time_step is None else int(time_step)   # see _ResolveTimeStep
        self.phase = float(phase)
        self._amplitude = None if amplitude is None else bool(amplitude)
        self.sample = max(int(sample), 1)
        self.time_sample = max(int(time_sample), 1)
        self.n_phase = max(int(n_phase), 1)
        self.animate = bool(animate)
        self.norm = ParseNorm(norm)
        self.unit = unit
        self._delta = src.delta_unit
        self._raw = {}           # read (normalised) data
        self._peaks = {}         # data key -> peak(s)
        self._step = None        # animation step, None: the static view
        self._anim = None
        self._anim_n = None

    @staticmethod
    def _Source(source, name, suffix=''):
        """ The FieldSource of `source` and the item name (default its stem). """
        from openEMS.Viewer.Source import AsSource
        src = AsSource(source)
        if name is None:
            if src.name is None:
                raise ValueError('a name is needed for a source from arrays')
            name = src.name + suffix
        return src, name

    def OnRemoved(self):
        # no open file while the item is not shown; read again, it reopens
        self.source.Close()

    def _ResolveTimeStep(self):
        """ Default time step: the one with the largest value (the first
        ones are often still empty, the last ones too after a pulse). Call
        once the subclass can read its data. """
        if self.time_step is not None:
            return
        n = self._NumTimes()
        if self.domain != 'time' or not n:
            self.time_step = n // 2
            return
        peaks = self._TimePeaks()
        self.time_step = max(peaks, key=peaks.get) if peaks else n // 2

    def _TimeSteps(self):
        """ The time indices of the animation, plus the static one. """
        n = self._NumTimes()
        steps = set(range(0, n, self.time_sample)) if n else set()
        if n and self.time_step is not None:
            steps.add(self._TimeIndex(None))
        return sorted(steps)

    def _TimePeaks(self):
        """ {time index: peak of the shown quantity} over :meth:`_TimeSteps`. """
        raise NotImplementedError

    # -- source access
    @property
    def amplitude(self):
        """ Frequency domain: the amplitude instead of the instantaneous value. """
        if self._amplitude is not None:
            return self._amplitude
        if self.AMPLITUDE_DEFAULT is not None:
            return self.AMPLITUDE_DEFAULT
        return self.quantity in ('magnitude', 'value')

    def _Kind(self):
        return self.source.kind

    def _CheckDomain(self, domain):
        if domain not in self.source.Domains():
            raise ValueError('domain "{}" not available, the dump offers {}'.format(domain, self.source.Domains()))

    def _SnapFreq(self, freq):
        """ A stored frequency nearest to `freq`, any one for a time-domain dump. """
        src = self.source
        if src.HasFD():
            return float(src.freqs[np.argmin(np.abs(src.freqs - float(freq)))])
        return float(freq)

    def _NumTimes(self):
        return len(self.source.times)

    def _TimeIndex(self, step):
        """ The time index shown: `time_step` in the static view, every
        `time_sample`-th one while animating. """
        n = self._NumTimes()
        if not n:
            return 0
        return self.time_step % n if step is None else (step * self.time_sample) % n

    def _Phase(self, step):
        return self.phase + 360.0 * (step or 0) / self.n_phase

    def _Norm(self, data):
        """ `data` divided by the normalisation factor (its magnitude for
        real time-domain data). """
        if self.norm == 1:
            return data
        if np.iscomplexobj(data):
            return data / self.norm
        return data / abs(self.norm)

    def _Unit(self):
        return 'dB' if self.scale == 'log' else (self.unit or self.source.unit)

    def _Invalidate(self):
        """ The data changed (quantity, domain, frequency, normalisation, ...). """
        self._raw = {}
        self._peaks = {}
        if self._anim is not None:
            self._anim.Invalidate()

    # -- texts
    def _WhenText(self, step=None, phase='φ='):
        """ The shown frequency (and phase) or time. """
        step = self._step if step is None else step
        if self.domain == 'frequency':
            txt = Quantity.FormatFreq(self.freq)
            if not self.amplitude:
                txt += ', {}{:.0f}°'.format(phase, self._Phase(step) % 360)
            return txt
        t = self._TimeIndex(step)
        return 't={}'.format(Quantity.FormatTime(self.source.times[t])) if self._NumTimes() else ''

    def _FixedWhenText(self):
        """ The shown frequency or time for texts that are not updated per
        frame (tree, info): without the phase or time while animated. """
        if self._anim is None:
            return self._WhenText()
        return Quantity.FormatFreq(self.freq) if self.domain == 'frequency' else ''

    # -- animation
    def _PhaseAnimated(self):
        """ Whether the frequency-domain view changes with the phase. """
        return not self.amplitude

    def _RegisterAnimation(self):
        v = self.GetViewer()
        if v is None:
            return
        n = None
        if self.animate:
            if self.domain == 'frequency':
                n = self.n_phase if self._PhaseAnimated() else None
            else:
                n = -(-self._NumTimes() // self.time_sample)
            if n is not None and n <= 1:
                n = None
        if n == self._anim_n and (n is None) == (self._anim is None):
            return
        if self._anim is not None:
            self._anim.Remove()
            self._anim = None
        self._anim_n = n
        if n is not None:
            first = v.GetAnimationSteps() == 0
            self._anim = self.AddAnimation(n)
            if first:
                self._SyncClock()

    def _SyncClock(self):
        """ Time domain: move the animation clock to the frame of
        `time_step`, so that the clock shows and plays from it. """
        v = self.GetViewer()
        if v is None or self._anim is None or self.domain != 'time':
            return
        n, N = self._anim.GetNumberOfSteps(), v.GetAnimationSteps()
        frame = int(round(self.time_step / self.time_sample)) % n
        # the first clock tick showing this frame
        v.SetAnimationStep(-(-frame * N // n), render=False)

    def _ShowStep(self, step):
        """ Show animation `step` (None: the static view); must not render. """
        raise NotImplementedError

    def OnAnimationStep(self, step):
        if self.domain == 'time':
            self.time_step = self._TimeIndex(step)
        self._ShowStep(step)

    def OnVisibilityChanged(self, visible):
        if visible:
            self._ShowStep(self._step)

    # -- settings
    def _FieldParameters(self):
        """ The settings of the shown field (section 'Field'). """
        src = self.source
        p = [dict(name='quantity', type='choice', value=self.quantity, choices=self._choices, section='Field')]
        if len(src.Domains()) > 1:
            p.append(dict(name='domain', type='choice', value=self.domain, choices=src.Domains(), section='Field'))
        if self.domain == 'frequency':
            if src.HasFD():
                p.append(dict(name='freq', label='frequency', type='choice', value=self.freq,
                              choices=[float(f) for f in src.freqs],
                              labels=[Quantity.FormatFreq(f) for f in src.freqs], section='Field'))
            else:
                p.append(dict(name='freq', label='frequency (Hz)', type='float', value=self.freq, min=0.0,
                              max=1e15, step=1e8, decimals=0, section='Field'))
            p.append(dict(name='phase', label='phase (°)', type='float', value=self.phase, min=-360.0,
                          max=360.0, step=15.0, default=0.0, section='Field'))
            if self.quantity != 'value':
                default = self.AMPLITUDE_DEFAULT
                if default is None:
                    default = self.quantity == 'magnitude'
                p.append(dict(name='amplitude', type='bool', value=self.amplitude, default=default,
                              section='Field'))
        else:
            p.append(dict(name='time_step', label='time step', type='int', value=self.time_step, min=0,
                          max=max(self._NumTimes() - 1, 0), section='Field'))
        return p

    def _DataParameters(self):
        """ Sampling and normalisation (section 'Data'). """
        norm = self.norm if isinstance(self.norm, float) else str(self.norm)
        return [dict(name='sample', label='spatial sampling', type='int', value=self.sample, min=1, max=100,
                     default=1, section='Data'),
                dict(name='norm', label='divide by', type='str', value=norm, default=1.0, section='Data'),
                dict(name='unit', label='unit', type='str', value=self.unit or '', default='', section='Data')]

    def _AnimationParameters(self):
        p = [dict(name='animate', type='bool', value=self.animate, default=True, section='Animation')]
        if self.domain == 'frequency':
            p.append(dict(name='n_phase', label='phase steps', type='int', value=self.n_phase, min=2, max=360,
                          default=36, section='Animation'))
        else:
            p.append(dict(name='time_sample', label='every n-th step', type='int', value=self.time_sample,
                          min=1, max=max(self._NumTimes(), 1), default=1, section='Animation'))
        return p

    def _SetField(self, name, value):
        """ Set a shared dump setting; False if `name` is none. """
        if name == 'quantity':
            self.quantity = Quantity.Normalize(value, self.source.kind, self._choices)
        elif name == 'domain':
            self._CheckDomain(value)
            if value == 'frequency' and self.freq is None:
                raise ValueError('set `freq` first')
            self.domain = value
        elif name == 'freq':
            self.freq = self._SnapFreq(value)
        elif name in ('time_step', 'sample', 'time_sample', 'n_phase'):
            setattr(self, name, max(int(value), 0 if name == 'time_step' else 1))
            if name == 'time_step':
                self._step = None
        elif name == 'phase':
            self.phase = float(value)
        elif name == 'amplitude':
            self._amplitude = bool(value)
        elif name == 'animate':
            self.animate = bool(value)
        elif name == 'norm':
            self.norm = ParseNorm(value)
        elif name == 'unit':
            self.unit = value or None
        else:
            return False
        return True

    def _SetOwn(self, name, value):
        """ Set a setting of the subclass; False if `name` is none. """
        return False

    def SetParameter(self, name, value):
        if not (self._SetField(name, value) or self._SetOwn(name, value) or self._SetCommon(name, value)):
            super().SetParameter(name, value)
        if name in self.DATA_SETTINGS:
            self._Invalidate()
        if self.GetViewer() is not None:
            self.Update()
            if name == 'time_step':
                self._SyncClock()
