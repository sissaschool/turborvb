#!/usr/bin/env bash

DESTINATION=../../src/m_common/fortran.c

python make_wrapper.py > $DESTINATION
