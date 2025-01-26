# Copyright (c) 2025.
# Authors: Ádám Fekete
# This program is distributed under the MIT License, see LICENSE.md.

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
