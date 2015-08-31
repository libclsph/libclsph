#!/bin/bash
find libclsph -type f | xargs clang-format-3.6 -style=file -i
