"""run a benchmark script, report its wall time and peak memory

/usr/bin/time is not in every image, and the runtime image has no build tools to add it.
"""
import resource, subprocess, sys, time

log = sys.argv[1]
start = time.time()
with open(log, 'w') as fh:
    subprocess.run([sys.executable, '-u'] + sys.argv[2:], stdout=fh, stderr=subprocess.STDOUT)
wall = time.time() - start
rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024   # Linux: KiB
print('  total: %.1f s' % wall)
print('  peak rss: %.0f MiB' % rss)
