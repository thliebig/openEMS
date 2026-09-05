.. _python_sar:

SAR
---

Specific Absorption Rate post-processing, see :ref:`concept_sar` for the
concept and the file formats.

Reading SAR Results
^^^^^^^^^^^^^^^^^^^

.. automodule:: openEMS.sar_utils

    .. autofunction:: readSAR

SAR Calculation
^^^^^^^^^^^^^^^

Computes local or mass-averaged SAR from a raw SAR dump
(``dump_type=29``). This is the same C++ implementation the ``sar_calc``
binary uses.

.. automodule:: openEMS.sar_calculation

    .. autoclass:: SAR_Calculation
        :members:
        :member-order: bysource
