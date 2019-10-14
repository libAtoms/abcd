class ABCDError(Exception):
    pass


class URLError(ABCDError):
    pass


class AuthenticationError(ABCDError):
    pass


class PropertyNotImplementedError(NotImplementedError):
    """Raised if a calculator does not implement the requested property."""


class TimeoutError(ABCDError):
    pass
