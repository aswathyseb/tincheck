"""
This module exists only to make module level execution also viable.

python -m tinrun

"""

from tinrun.main import run


def main():
    """
    A simple wrapper over the task_selector
    """
    run()


if __name__ == '__main__':
    main()
