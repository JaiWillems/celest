"""Define decorators used throughout Celest."""


def set_module(module):
    """Override the module of a class or function."""

    def decorator(func):

        if module is not None:
            func.__module__ = module

        return func

    return decorator
