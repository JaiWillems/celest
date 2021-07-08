"""Define decorators used throughout Celest."""


def set_module(module):
    """Specify a class or functions module."""

    def decorator(func):

        if module is not None:
            func.__module__ = module

        return func

    return decorator