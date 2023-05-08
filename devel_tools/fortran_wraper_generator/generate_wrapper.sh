#!/usr/bin/env bash

DESTINATION=../../src/m_common/fortran.c

python3 make_wrapper.py > $DESTINATION
