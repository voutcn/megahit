#!/usr/bin/env python3

import random


def main():
    random.seed(10086)
    ACGT = "ACGT"
    for i in range(1000):
        print('>{}'.format(i + 1))
        print(''.join(ACGT[random.randint(0, 3)] for _ in range(300)))


if __name__ == '__main__':
    main()
