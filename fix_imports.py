with open('data_processing/tests.py', 'r') as f:
    lines = f.readlines()

new_lines = []
imports = []
code = []

in_import = False
for line in lines:
    if line.startswith('from data_processing'):
        imports.append(line)
        in_import = True
    elif in_import:
        imports.append(line)
        if line.strip() == ')':
            in_import = False
    elif line.startswith('import unittest') or line.startswith('import pandas') or line.startswith('from unittest.mock') or line.startswith('import logging'):
        imports.append(line)
    elif line.startswith('logging.basicConfig'):
        code.append(line)
    else:
        code.append(line)

with open('data_processing/tests.py', 'w') as f:
    f.writelines(imports + code)
