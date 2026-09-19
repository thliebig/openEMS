# test-only: run every openEMS simulation with the engine given in OPENEMS_TEST_ENGINE
import os, sys
_engine = os.environ.get('OPENEMS_TEST_ENGINE')
if _engine:
    import openEMS
    _mod = sys.modules['openEMS.openEMS']
    class _EngineOverride(_mod.openEMS):
        def Run(self, sim_path, *args, **kw):
            kw.setdefault('engine', _engine)
            return super().Run(sim_path, *args, **kw)
    openEMS.openEMS = _EngineOverride
    _mod.openEMS = _EngineOverride
