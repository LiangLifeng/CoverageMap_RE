"""Exclusive process CPU/wall accounting, including nested stages and exceptions."""
from contextlib import contextmanager
import resource
import time


class StageTiming:
    def __init__(self):
        self.values = {}
        self.active = 'other'
        self.previous = self.clock()

    @staticmethod
    def clock():
        usage = resource.getrusage(resource.RUSAGE_SELF)
        return time.perf_counter(), usage.ru_utime, usage.ru_stime

    def charge(self):
        now = self.clock()
        row = self.values.setdefault(self.active, dict(wall_seconds=0., user_cpu_seconds=0.,
                                                       system_cpu_seconds=0., calls=0))
        for key, old, new in zip(('wall_seconds', 'user_cpu_seconds', 'system_cpu_seconds'), self.previous, now):
            row[key] += max(0., new - old)
        self.previous = now

    @contextmanager
    def measure(self, name):
        self.charge()
        parent, self.active = self.active, name
        self.values.setdefault(name, dict(wall_seconds=0., user_cpu_seconds=0., system_cpu_seconds=0., calls=0))['calls'] += 1
        try:
            yield
        finally:
            self.charge()
            self.active = parent

    def snapshot(self):
        self.charge()
        return dict(schema='exclusive-process-stages-1', stages={k: dict(v) for k, v in sorted(self.values.items())},
                    scope='process CPU includes its threads; exclusive elapsed measurement windows; thread CPU may overlap output/callback work; excludes snapshot publication and CLI startup')
