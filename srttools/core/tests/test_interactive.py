from __future__ import (absolute_import, division,
                        print_function)
import numpy as np
from ..interactive_filter import ImageSelector, DataSelector
import warnings
import pytest


np.random.seed(1241347)


class TestImageSelector(object):
    @classmethod
    def setup_class(klass):
        pass

    def test_interactive_valid_data(self):
        import matplotlib.pyplot as plt
        data = np.zeros((100, 100))
        ax = plt.subplot()
        a = ImageSelector(data, ax, test=True)
        fake_event = type('event', (), {})()
        fake_event.key = 'q'
        fake_event.xdata, fake_event.ydata = (130, 30)

        retval = a.on_key(fake_event)
        assert retval == (130, 30, 'q')

    def test_interactive_invalid_data(self):
        import matplotlib.pyplot as plt
        data = np.zeros((100, 100))
        ax = plt.subplot()
        a = ImageSelector(data, ax, test=True)
        fake_event = type('event', (), {})()
        fake_event.key = 'b'
        fake_event.xdata, fake_event.ydata = (None, 30)

        retval = a.on_key(fake_event)
        assert retval is None

    def test_interactive_fun(self):
        import matplotlib.pyplot as plt
        data = np.zeros((100, 100))
        ax = plt.subplot()
        def fun(x, y, key):
            warnings.warn("It is working", UserWarning)
        a = ImageSelector(data, ax, test=True, fun=fun)
        fake_event = type('event', (), {})()
        fake_event.key = 'b'
        fake_event.xdata, fake_event.ydata = (130, 30)

        with pytest.warns(UserWarning) as record:
            retval = a.on_key(fake_event)
        assert "It is working" in record[0].message.args[0]
        assert retval == (130, 30, 'b')
