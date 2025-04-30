"""
This is a minimalistic performance timer.
"""
import time

__author__ = "Gregory Halverson"

DEFAULT_FORMAT = "0.2f"

class Timer(object):
    """
    This is a minimalistic performance timer.
    """

    def __init__(self):
        self._start_time = None
        self._end_time = None
        self.start()

    def __enter__(self, *args, **kwargs):
        self.start()
        return self

    def __exit__(self, *args, **kwargs):
        self.end()

    def __repr__(self):
        # print("Timer.__repr__")
        return self.__format__(format_string=DEFAULT_FORMAT)

    def __str__(self):
        # print("Timer.__str__")
        return self.__repr__()

    def __format__(self, format_string=DEFAULT_FORMAT):
        if format_string is None or format_string == "":
            format_string = DEFAULT_FORMAT

        return format(self.duration, format_string)

    @property
    def now(self):
        # return datetime.now()
        return time.perf_counter()

    def start(self):
        self._start_time = self.now

        return self.start_time

    @property
    def start_time(self):
        return self._start_time

    def end(self):
        self._end_time = self.now

        return self.end_time

    @property
    def end_time(self):
        return self._end_time

    @property
    def duration(self):
        if self.start_time is None:
            raise Exception("timer never started")

        if self.end_time is None:
            end_time = self.now
        else:
            end_time = self.end_time

        duration = end_time - self.start_time

        return duration

