import numpy as np
from collections import Counter

from abcd import ABCD


if __name__ == '__main__':
    abcd = ABCD(url='mongodb://localhost:27017')
    print(abcd)
    abcd.print_info()
