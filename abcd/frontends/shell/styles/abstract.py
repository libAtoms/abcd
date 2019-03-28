from abc import ABCMeta, abstractmethod


class Style(metaclass=ABCMeta):
    # @abstractmethod
    # def __init__(self, *args, **kwargs):
    #     pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    @abstractmethod
    def title(self, string):
        pass

    @abstractmethod
    def h1(self, string):
        pass

    @abstractmethod
    def h2(self, string):
        pass

    @abstractmethod
    def hist(self, data, **kwargs):
        pass

    @abstractmethod
    def hist_str(self, data, max_width=80):
        pass

    @abstractmethod
    def hist_float(self, data, bins=10):
        pass

    @staticmethod
    @abstractmethod
    def print(*args, **kwargs):
        pass

    @staticmethod
    @abstractmethod
    def table(*args, **kwargs):
        pass

    @staticmethod
    @abstractmethod
    def newline():
        pass
