"""
A Python Singleton mixin class that makes use of some of the ideas
found at http://c2.com/cgi/wiki?PythonSingleton. Just inherit
from it and you have a singleton. No code is required in
subclasses to create singleton behavior -- inheritance from 
Singleton is all that is needed.

Singleton creation is threadsafe.

USAGE:

Just inherit from Singleton. If you need a constructor, include
an __init__() method in your class as you usually would. However,
if your class is S, you instantiate the singleton using S.getInstance() 
instead of S(). Repeated calls to S.getInstance() return the 
originally-created instance.

For example:

class S(Singleton):

    def __init__(self, a, b=1):
        pass
        
S1 = S.getInstance(1, b=3)


Most of the time, that's all you need to know. However, there are some
other useful behaviors. Read on for a full description:

1) Getting the singleton:

    S.getInstance() 
    
returns the instance of S. If none exists, it is created. 

2) The usual idiom to construct an instance by calling the class, i.e.

    S()
    
is disabled for the sake of clarity. 

For one thing, the S() syntax means instantiation, but getInstance()
usually does not cause instantiation. So the S() syntax would
be misleading.

Because of that, if S() were allowed, a programmer who didn't 
happen to notice the inheritance from Singleton (or who
wasn't fully aware of what a Singleton pattern
does) might think he was creating a new instance, 
which could lead to very unexpected behavior.

So, overall, it is felt that it is better to make things clearer
by requiring the call of a class method that is defined in
Singleton. An attempt to instantiate via S() will result 
in a SingletonException being raised.

3) Use __S.__init__() for instantiation processing,
since S.getInstance() runs S.__init__(), passing it the args it has received. 

If no data needs to be passed in at instantiation time, you don't need S.__init__().

4) If S.__init__(.) requires parameters, include them ONLY in the
first call to S.getInstance(). If subsequent calls have arguments,
a SingletonException is raised by default.

If you find it more convenient for subsequent calls to be allowed to
have arguments, but for those argumentsto be ignored, just include 
'ignoreSubsequent = True' in your class definition, i.e.:

  class S(Singleton):
  
      ignoreSubsequent = True

      def __init__(self, a, b=1):
          pass

5) For testing, it is sometimes convenient for all existing singleton
instances to be forgotten, so that new instantiations can occur. For that
reason, a forgetAllSingletons() function is included. Just call

  forgetAllSingletons()
  
and it is as if no earlier instantiations have occurred.

6) As an implementation detail, classes that inherit 
from Singleton may not have their own __new__
methods. To make sure this requirement is followed, 
an exception is raised if a Singleton subclass includ
es __new__. This happens at subclass instantiation
time (by means of the MetaSingleton metaclass.


By Gary Robinson, grobinson@flyfi.com. No rights reserved -- 
placed in the public domain -- which is only reasonable considering
how much it owes to other people's code and ideas which are in the
public domain. The idea of using a metaclass came from 
a  comment on Gary's blog (see 
http://www.garyrobinson.net/2004/03/python_singleto.html#comments). 
Other improvements came from comments and email from other
people who saw it online. (See the blog post and comments
for further credits.)

Not guaranteed to be fit for any particular purpose. Use at your
own risk. 
"""

import threading


class SingletonException(Exception):
    pass


_SINGLETONS = set()
_LOCK = threading.RLock()

# Ensure only one instance of each Singleton class is created.  This is not bound to the
# individual Singleton class since we need to ensure that there is only one mutex for each
# Singleton class, which would require having a lock when setting up the Singleton class,
# which is what this is anyway. So, when any # Singleton is created, we lock this lock and
# then we don't need to lock it again for that class.
_CREATE_LOCK = threading.RLock()


def _create_singleton_instance(cls, args, kwargs):
    _CREATE_LOCK.acquire()
    try:
        if cls._is_instantiated():  # some other thread got here first
            return 
        
        instance = cls.__new__(cls)
        try:
            instance.__init__(*args, **kwargs)
        except TypeError as e:
            if str(e).find('__init__() takes') != -1:
                msg = 'If the singleton requires __init__ args, supply them on first call to get_instance().'
                raise SingletonException(msg)
            else:
                raise
        cls.c_instance = instance
        _add_singleton(cls)
    finally:
        _CREATE_LOCK.release()


def _add_singleton(cls):
    _LOCK.acquire()
    try:
        assert cls not in _SINGLETONS
        _SINGLETONS.add(cls)
    finally:
        _LOCK.release()


def _remove_singleton(cls):
    _LOCK.acquire()
    try:
        if cls in _SINGLETONS:
            _SINGLETONS.remove(cls)
    finally:
        _LOCK.release()


def forget_all_singletons():
    """

    This is useful in tests, since it is hard to know which singletons need to be cleared to make a test work.

    """
    _LOCK.acquire()
    try:
        for cls in _SINGLETONS.copy():
            cls._forget_class_instance_reference()

        # Might have created some Singletons in the process of tearing down.
        # Try one more time - there should be a limit to this.
        number_of_singletons = len(_SINGLETONS)
        if len(_SINGLETONS) > 0:
            for cls in _SINGLETONS.copy():
                cls._forget_class_instance_reference()
                number_of_singletons -= 1
                assert number_of_singletons == len(_SINGLETONS), 'Added a singleton while destroying ' + str(cls)
        assert len(_SINGLETONS) == 0, _SINGLETONS
    finally:
        _LOCK.release()


class MetaSingleton(type):
    def __new__(mcs, name, bases: tuple, attributes: dict):
        if '__new__' in attributes:
            raise SingletonException('Can not override __new__ in a Singleton')
        return super(MetaSingleton, mcs).__new__(mcs, name, bases, attributes)
        
    def __call__(cls, *args, **kwargs):
        raise SingletonException('Singletons may only be instantiated through get_instance()')


class Singleton(object, metaclass=MetaSingleton):
    def get_instance(cls, *args, **kwargs):
        """
        Call this to instantiate an instance or retrieve the existing instance.
        If the singleton requires args to be instantiated, include them the first
        time you call getInstance.        
        """
        if cls._is_instantiated():
            if (args or kwargs) and not hasattr(cls, 'ignore_subsequent'):
                raise SingletonException('Singleton already instantiated, but get_instance() called with args.')
        else:
            _create_singleton_instance(cls, args, kwargs)
            
        return cls.c_instance
    get_instance = classmethod(get_instance)
    
    def _is_instantiated(cls):
        # Don't use hasattr(cls, 'cInstance'), because that screws things up if there is a singleton that
        # extends another singleton.  hasattr looks in the base class if it doesn't find in subclass.
        return 'c_instance' in cls.__dict__
    _is_instantiated = classmethod(_is_instantiated)

    # This can be handy for public use also
    is_instantiated = _is_instantiated

    def _forget_class_instance_reference(cls):
        """
        This is designed for convenience in testing -- sometimes you 
        want to get rid of a singleton during test code to see what
        happens when you call getInstance() under a new situation.
        
        To really delete the object, all external references to it
        also need to be deleted.
        """
        try:
            if hasattr(cls.c_instance, '_prepare_to_forget_singleton'):
                # tell instance to release anything it might be holding onto.
                cls.c_instance._prepare_to_forget_singleton()
            del cls.c_instance
            _remove_singleton(cls)
        except AttributeError:
            # run up the chain of base classes until we find the one that has the instance
            # and then delete it there
            for baseClass in cls.__bases__: 
                if issubclass(baseClass, Singleton):
                    baseClass._forget_class_instance_reference()
    _forget_class_instance_reference = classmethod(_forget_class_instance_reference)
