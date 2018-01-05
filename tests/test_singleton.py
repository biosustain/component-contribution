import threading

import pytest
import time

from component_contribution.core.singleton import Singleton, SingletonException, MetaSingleton, \
    forget_all_singletons


def test_singleton_returns_the_same_object():
    """
    Demonstrates normal use -- just call get_instance and it returns a singleton instance
    """

    class A(Singleton):
        def __init__(self):
            super(A, self).__init__()

    a1 = A.get_instance()
    a2 = A.get_instance()
    assert id(a1) == id(a2)

def test_instantiate_using_arguments():
    """
    If the singleton needs args to construct, include them in the first
    call to get instances.
    """

    class B(Singleton):
        def __init__(self, arg1, arg2):
            super(B, self).__init__()
            self.arg1 = arg1
            self.arg2 = arg2

    b1 = B.get_instance('arg1 value', 'arg2 value')
    b2 = B.get_instance()
    assert b1.arg1 == 'arg1 value'
    assert b1.arg2 == 'arg2 value'
    assert id(b1) == id(b2)


def test_instantiate_using_keyword_arguments():

    class B(Singleton):
        def __init__(self, arg1=5):
            super(B, self).__init__()
            self.arg1 = arg1

    b1 = B.get_instance('arg1 value')
    b2 = B.get_instance()
    assert b1.arg1 == 'arg1 value'
    assert id(b1) == id(b2)


def test_instantiate_without_required_arguments():

    class B(Singleton):
        def __init__(self, arg1, arg2):
            super(B, self).__init__()
            self.arg1 = arg1
            self.arg2 = arg2

    with pytest.raises(TypeError):
        B.get_instance()


def test_pass_type_error():
    """
    Make sure the test for capturing missing args doesn't interfere with a normal TypeError.
    """

    class B(Singleton):
        def __init__(self, arg1, arg2):
            super(B, self).__init__()
            self.arg1 = arg1
            self.arg2 = arg2
            raise TypeError('some type error')

    with pytest.raises(TypeError):
        B.get_instance(1, 2)


def test_instantiate_without_get_instance():
    """
    Demonstrates that singletons can ONLY be instantiated through
    getInstance, as long as they call Singleton.__init__ during construction.

    If this check is not required, you don't need to call Singleton.__init__().
    """

    class A(Singleton):
        def __init__(self):
            super(A, self).__init__()

    with pytest.raises(SingletonException):
        A()


def test_dont_allow_new():

    def instantiated_illegal_class():
        class A(Singleton):
            def __init__(self):
                super(A, self).__init__()

            def __new__(mtc, name, bases, attributes):
                return super(MetaSingleton, mtc).__new__(mtc, name, bases, attributes)

    with pytest.raises(SingletonException):
        instantiated_illegal_class()


def test_dont_allow_args_after_construction():
    class B(Singleton):
        def __init__(self, arg1, arg2):
            super(B, self).__init__()
            self.arg1 = arg1
            self.arg2 = arg2

    B.get_instance('arg1 value', 'arg2 value')
    with pytest.raises(SingletonException):
        B('arg1 value', 'arg2 value')


def test_forget_class_instance_reference():
    class A(Singleton):
        def __init__(self):
            super(A, self).__init__()

    class B(A):
        def __init__(self):
            super(B, self).__init__()

    # check that changing the class after forgetting the instance produces
    # an instance of the new class
    a = A.get_instance()
    assert a.__class__.__name__ == 'A'
    A._forget_class_instance_reference()
    b = B.get_instance()
    assert b.__class__.__name__ == 'B'

    # check that invoking the 'forget' on a subclass still deletes the instance
    B._forget_class_instance_reference()
    a = A.get_instance()
    B._forget_class_instance_reference()
    b = B.get_instance()
    assert b.__class__.__name__ == 'B'


def test_forget_all_singletons():
    # Should work if there are no singletons
    forget_all_singletons()

    class A(Singleton):
        init_count = 0

        def __init__(self):
            super(A, self).__init__()
            A.init_count += 1

    A.get_instance()
    assert A.init_count == 1

    A.get_instance()
    assert A.init_count == 1

    forget_all_singletons()
    A.get_instance()
    assert A.init_count == 2


def test_threaded():
    """
    Check that only one Singleton is created even if multiple
    threads try at the same time.  If fails, would see assert in _add_singleton
    """

    class TestSingleton(Singleton):
        def __init__(self):
            super(TestSingleton, self).__init__()

    class TestSingletonThread(threading.Thread):
        def __init__(self, target_time):
            super(TestSingletonThread, self).__init__()
            self._target_time = target_time
            self._exception = None

        def run(self):
            try:
                sleep_time = self._target_time - time.time()
                if sleep_time > 0:
                    time.sleep(sleep_time)
                TestSingleton.get_instance()
            except Exception as e:
                self._exception = e

    target_time = time.time() + 0.1
    threads = []
    for _ in range(100):
        t = TestSingletonThread(target_time)
        t.start()
        threads.append(t)
    exception = None
    for t in threads:
        t.join()
        if t._exception and not exception:
            exception = t._exception
    assert exception is None


def test_no_init():
    """
    Demonstrates use with a class not defining __init__
    """

    class A(Singleton):
        pass

    A.get_instance()  # Make sure no exception is raised


def test_multiple_get_instances_with_arguments():

    class A(Singleton):
        ignore_subsequent = True

        def __init__(self, a, b=1):
            pass

    A.get_instance(1)
    A.get_instance(2)  # ignores the second call because of ignoreSubsequent

    class B(Singleton):
        def __init__(self, a, b=1):
            pass

    B.get_instance(1)
    with pytest.raises(SingletonException):
        B.get_instance(2)  # No ignoreSubsequent included

    class C(Singleton):
        def __init__(self, a=1):
            pass

    C.get_instance(a=1)
    with pytest.raises(SingletonException):
        C.get_instance(a=2)  # No ignoreSubsequent included


def test_inheritance():
    """
    It's sometimes said that you can't subclass a singleton (see, for instance,
    http://steve.yegge.googlepages.com/singleton-considered-stupid point e). This
    test shows that at least rudimentary subclassing works fine for us.
    """

    class A(Singleton):
        def __init__(self):
            self._x = None
            self._z = None

        @property
        def x(self):
            return self._x

        @x.setter
        def x(self, x):
            self._x = x

        @property
        def z(self):
            return self._z

        @z.setter
        def z(self, z):
            raise NotImplementedError

    class B(A):
        def __init__(self):
            super(A, self).__init__()
            self._y = None

        @A.x.setter
        def x(self, x):
            self._x = -x

        @property
        def y(self):
            return self._y

        @y.setter
        def y(self, y):
            self._y = y

    a = A.get_instance()
    a.x = 5
    b = B.get_instance()
    b.x = 5
    b.y = 50
    assert (a.x, b.x, b.y) == (5, -5, 50)
    with pytest.raises(AttributeError):
        a.y
    with pytest.raises(NotImplementedError):
        b.z = 500
