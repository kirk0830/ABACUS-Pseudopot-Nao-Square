"""this file is to run all unittest within the apns package"""
import os

def run_all_tests():
    """run all unittests in files within current and all subfolders.
    Present all unittests are written behind the implementation of
    the module in the same file.
    """
    # get the current directory
    cwd = os.path.dirname(__file__)
    # simply run with `python3 xxx.py` is enough for running unitests
    # programed in the same file as the module

    for root, dirs, files in os.walk(cwd):
        for file in files:
            # if there is `unittest` imported, then run the unittest, otherwise
            # print warning that file is not programed with unittest
            if file.endswith('.py'):
                with open(os.path.join(root, file), 'r') as f:
                    content = f.read()
                    if 'import unittest' in content:
                        print('Running unittest in', file)
                        os.system('python3 ' + os.path.join(root, file))
                    else:
                        print('Warning: no unittest in', file)

if __name__ == '__main__':
    """run all unittests in the apns package. This will be equipped in Github
    automatic testing."""
    run_all_tests()