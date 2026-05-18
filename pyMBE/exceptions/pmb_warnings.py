import warnings
from functools import wraps

def deprecated(new_function):
    """Wrapper function to display a warning deprecation message to all functions that are going to be removed in a future
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            warnings.warn(
                f"{func.__name__} is deprecated and it will be removed in the future. Please use"+ new_function +"instead.",
                DeprecationWarning,
                stacklevel=2
            )
            return func(*args, **kwargs)
        return wrapper
    return decorator
